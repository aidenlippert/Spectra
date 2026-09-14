"""Bounded full-spectrum soft-min annealing; untrusted proposal only."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,sys,time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from joint_profile_numeric import build,profiles
from experiments.marginal_charged_projectors import charged_vectors
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value
from signed_charge_numeric import prepare
from experiments.marginal_hopping_telescope import actions,projected_matrix
from experiments.marginal_spin_telescope import LABELS,actions as spin_actions
from experiments.marginal_spectator_hopping import LABELS as SPECTATORS,actions as spectator_actions
from experiments.marginal_pair_transfer import LABELS as PAIRS,actions as pair_actions
from experiments.marginal_two_spectator_hopping import LABELS as TWO,actions as two_actions
from experiments.marginal_three_spectator_hopping import LABELS as THREE,actions as three_actions
from experiments.marginal_coherent_projector_telescope import LABELS as COHERENT,actions as coherent_actions
BASE=ROOT/'results/marginal_graded_hubbard8'



def main():
    parser=argparse.ArgumentParser();parser.add_argument('seed',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();c=json.loads(args.seed.read_text())
    indices=[1,2,4,5,6,7,31,11,20]
    candidates_path=BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json';candidates=json.loads(candidates_path.read_text())['candidates']
    shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in indices]
    x,physical,mats,*_=prepare(c,shapes)
    action=[actions()]+[spin_actions({label:1}) for label in LABELS]+[spectator_actions({label:1}) for label in SPECTATORS]+[pair_actions({label:1}) for label in PAIRS]+[two_actions({label:1}) for label in TWO]+[three_actions({label:1}) for label in THREE]+[coherent_actions({label:1}) for label in COHERENT];hop={}
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);hop[key]=np.array([projected_matrix(cols,item) for item in action],float)/scale[:,None]/scale[None,:]
    theta_h=F(c['projector_sum_ceiling'])/c['windows'];theta_j=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
    diagonal={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};taus=[]
    for shape in shapes:
        ratios={diagonal.get(s,F(0))/v for s,v in shape.items()}
        if len(ratios)!=1:raise ValueError('Seed sparse correction is outside the fixed span')
        taus.append(float(ratios.pop()))
    y0=np.array(list(map(float,x[2:4]))+[float(F(c['penalty'])),float(F(c['joint']['penalty']))]+taus+[float(F(c.get('signed_charge_telescope',{}).get(key,0))) for key in PATTERNS]+[float(F(c.get('hopping_telescope',0)))]+[float(F(c.get('spin_telescope',{}).get(label,0))) for label in LABELS]+[float(F(c.get('spectator_hopping',{}).get(label,0))) for label in SPECTATORS]+[float(F(c.get('pair_transfer',{}).get(label,0))) for label in PAIRS]+[float(F(c.get('two_spectator_hopping',{}).get(label,0))) for label in TWO]+[float(F(c.get('three_spectator_hopping',{}).get(label,0))) for label in THREE]+[float(F(c.get('coherent_projector',{}).get(label,0))) for label in COHERENT]+[float(F(c['penalized_lower']))])
    assert len(y0)==139
    cached=None;cached_tau=None;value=None;gradient=None;best=None;calls=0;tau=1e-5
    stages=[];bound_checks=0;latest_spectra=None
    def evaluate(y):
        nonlocal cached,cached_tau,value,gradient,best,calls,bound_checks,latest_spectra
        if cached is not None and cached_tau==tau and np.array_equal(y,cached):return value,gradient
        if calls>=499:raise RuntimeError('Optimizer499 evaluation budget reached; one reserved for atom export')
        calls+=1;cached=y.copy();cached_tau=tau;spectra=[];minimum=float('inf')
        for key,a,ds,ph,q,ts,*_ in mats:
            matrix=a+np.einsum('i,ijk->jk',y[:2]-np.array(x[2:4],float),ds)+y[2]*ph+y[3]*q+np.diag(y[4:65]@ts)+np.einsum('i,ijk->jk',y[65:138],hop[key])
            ev,vectors=eigh(matrix);minimum=min(minimum,float(ev[0]));spectra.append((key,ev,vectors,ds,ph,q,ts))
        latest_spectra=spectra
        partition=sum(float(np.exp(-(ev-minimum)/tau).sum()) for key,ev,v,ds,ph,q,ts in spectra)
        smooth=minimum-tau*math.log(partition)
        if not minimum-tau*math.log(4096)-1e-12<=smooth<=minimum+1e-12:raise ValueError('Soft-min bound failed')
        bound_checks+=1
        expectation=np.zeros(138)
        for key,ev,vectors,ds,ph,q,ts in spectra:
            weights=np.exp(-(ev-minimum)/tau)/partition
            if not np.any(weights):continue
            rho=(vectors*weights)@vectors.T
            expectation+=np.r_[[np.sum(rho*d) for d in ds],np.sum(rho*ph),np.sum(rho*q),ts@np.diag(rho),np.einsum('ij,kij->k',rho,hop[key])]
        gradient=-expectation/5;gradient[2]+=float(theta_h)/5;gradient[3]+=float(theta_j)/5
        value=(y[2]*float(theta_h)+y[3]*float(theta_j)-smooth)/5
        score=(minimum-y[2]*float(theta_h)-y[3]*float(theta_j))/5
        if min(y[:4])>=0 and y[0]+y[1]<=2.5+1e-12 and max(abs(y[13:65]))<=8+1e-9 and max(abs(np.r_[y[4:13],y[65:138]]))<=2+1e-9 and (best is None or score>best[0]):best=(score,np.r_[y,minimum])
        return value,gradient
    point=y0[:-1].copy();initial_value,initial_gradient=evaluate(point)
    rng=np.random.default_rng(20260912);gradient_errors=[]
    for _ in range(2):
        direction=rng.normal(size=138);direction/=np.linalg.norm(direction);step=1e-8
        plus=evaluate(point+step*direction)[0];minus=evaluate(point-step*direction)[0]
        error=abs((plus-minus)/(2*step)-initial_gradient@direction);gradient_errors.append(float(error))
        if error>2e-7:raise ValueError('Smoothed gradient finite-difference check failed')
    curvature_checks=[]
    def coordinate_scale(y):
        objective,g=evaluate(y)
        minimum=min(float(ev[0]) for key,ev,v,ds,ph,q,ts in latest_spectra)
        partition=sum(float(np.exp(-(ev-minimum)/tau).sum()) for key,ev,v,ds,ph,q,ts in latest_spectra)
        curvature=np.zeros((138,138))
        for key,ev,vectors,ds,ph,q,ts in latest_spectra:
            weights=np.exp(-(ev-minimum)/tau)/partition
            for i in np.flatnonzero(weights>1e-14):
                v=vectors[:,i]
                products=np.column_stack([np.column_stack([d@v for d in ds]),ph@v,q@v,ts.T*v[:,None],np.column_stack([a@v for a in hop[key]])])
                projected=vectors.T@products
                curvature+=weights[i]*np.outer(projected[i],projected[i])/tau
                gaps=ev[i+1:]-ev[i]
                coefficients=np.empty(len(gaps));near=gaps<tau*1e-8
                coefficients[near]=weights[i]/tau
                coefficients[~near]=weights[i]*(-np.expm1(-gaps[~near]/tau))/gaps[~near]
                curvature+=2*projected[i+1:].T@(coefficients[:,None]*projected[i+1:])
        mean=-5*g;mean[2]+=float(theta_h);mean[3]+=float(theta_j)
        hessian=(curvature-np.outer(mean,mean)/tau)/5
        diagonal=np.diag(hessian)
        # Check two coordinate curvatures against a derivative difference.
        errors=[]
        for j in (60,70):
            step=1e-9;direction=np.zeros(138);direction[j]=step
            numeric=(evaluate(y+direction)[1][j]-evaluate(y-direction)[1][j])/(2*step)
            error=abs(numeric-diagonal[j]);errors.append(float(error))
            if error>1e-2*(1+abs(diagonal[j])):raise ValueError('Spectral curvature check failed')
        direction=np.random.default_rng(20260912).normal(size=138);direction/=np.linalg.norm(direction);step=1e-9
        numeric=(evaluate(y+step*direction)[1]-evaluate(y-step*direction)[1])/(2*step)
        expected=hessian@direction;direction_error=float(max(abs(numeric-expected)))
        if direction_error>.02*(1+max(abs(expected))):raise ValueError('Full spectral Hessian check failed')
        eigenvalues,eigenvectors=np.linalg.eigh((hessian+hessian.T)/2)
        scales=np.clip(1/np.sqrt(np.maximum(eigenvalues,1e-10)),1e-3,100.)
        transform=(eigenvectors*scales)@eigenvectors.T
        curvature_checks.append({'temperature':tau,'minimum_curvature':float(min(diagonal)),'maximum_curvature':float(max(diagonal)),'scale_min':float(min(scales)),'scale_max':float(max(scales)),'finite_difference_errors':errors,'directional_hessian_error':direction_error,'hessian_min_eigenvalue':float(min(eigenvalues)),'hessian_max_eigenvalue':float(max(eigenvalues))})
        return transform
    outcome=None
    for tau in (1e-6,1e-7):
        before=calls
        try:
            origin=point.copy();transform=coordinate_scale(point)
            bounds=np.array([(0,2.5)]*2+[(0,3)]*2+[(-2,2)]*9+[(-8,8)]*52+[(-2,2)]*73)
            def constraints(z):
                y=origin+transform@z
                return np.r_[y-bounds[:,0],bounds[:,1]-y,2.5-y[0]-y[1]]
            constraint_jac=np.vstack([transform,-transform,-transform[:2].sum(axis=0)])
            outcome=minimize(lambda z:evaluate(origin+transform@z)[0],np.zeros(138),jac=lambda z:transform.T@evaluate(origin+transform@z)[1],method='SLSQP',constraints=[{'type':'ineq','fun':constraints,'jac':lambda z:constraint_jac}],options={'maxiter':120,'ftol':1e-13})
            point=origin+transform@outcome.x
            stages.append({'temperature':tau,'success':bool(outcome.success),'message':str(outcome.message),'evaluations':calls-before,'best_periodic_lower':best[0]})
        except RuntimeError as error:
            if 'budget' not in str(error):raise
            stages.append({'temperature':tau,'success':False,'message':str(error),'evaluations':calls-before,'best_periodic_lower':best[0]});outcome=None;break
        print(json.dumps(stages[-1]),flush=True)
    # One reserved spectral evaluation exports untrusted physical directions.
    atom_point=best[1][:-1];atom_spectra=[];atom_minimum=float('inf');calls+=1
    if calls>500:raise RuntimeError('Total500 full-spectrum budget exceeded')
    for key,a,ds,ph,q,ts,*_ in mats:
        matrix=a+np.einsum('i,ijk->jk',atom_point[:2]-np.array(x[2:4],float),ds)+atom_point[2]*ph+atom_point[3]*q+np.diag(atom_point[4:65]@ts)+np.einsum('i,ijk->jk',atom_point[65:138],hop[key])
        ev,vectors=eigh(matrix);atom_minimum=min(atom_minimum,float(ev[0]));atom_spectra.append((key,ev,vectors))
    z=sum(float(np.exp(-(ev-atom_minimum)/tau).sum()) for key,ev,v in atom_spectra)
    atoms=[];omitted=0.
    for key,ev,vectors in atom_spectra:
        weights=np.exp(-(ev-atom_minimum)/tau)/z
        for i,weight in enumerate(weights):
            if weight>1e-12:atoms.append({'sector':key,'dual_weight':float(weight)/5,'coordinates':vectors[:,i].tolist()})
            else:omitted+=float(weight)
    (args.output/'dual_atoms.json').write_text(json.dumps({'accepted':False,'atoms':atoms,'temperature':tau,'omitted_weight':omitted,'evaluation_point':atom_point.tolist(),'scope':'Untrusted Gibbs-state directions at the best unsmoothed point; no exact stationarity or family ceiling claim.'})+'\n')
    (args.output/'thermal_history.json').write_text(json.dumps({'stages':stages,'gradient_errors':gradient_errors,'soft_min_bound_checks':bound_checks,'best_unrounded_point':best[1].tolist(),'last_optimizer_point':point.tolist(),'curvature_checks':curvature_checks},indent=2)+'\n')
    y=[F(round(float(value)*10**6),10**6) for value in best[1][:-1]]
    point=x.copy();point[2:4]=y[:2];fresh=physical(point);signed=dict(zip(PATTERNS,y[13:65]));sparse={s:sum(t*shape.get(s,0) for t,shape in zip(y[4:13],shapes)) for s in set().union(*shapes)}
    fresh_action=[actions()]+[spin_actions({label:1}) for label in LABELS]+[spectator_actions({label:1}) for label in SPECTATORS]+[pair_actions({label:1}) for label in PAIRS]+[two_actions({label:1}) for label in TWO]+[three_actions({label:1}) for label in THREE]+[coherent_actions({label:1}) for label in COHERENT]
    minimum=(float('inf'),None);difference=0.
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        k=fresh[key][1];scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);direct=np.array(k,float)/scale[:,None]/scale[None,:]+np.einsum('i,ijk->jk',np.array(y[65:138],float),np.array([projected_matrix(cols,item) for item in fresh_action],float))/scale[:,None]/scale[None,:]
        correction=np.array([float(sparse.get(next(iter(col))&1023,0)-sparse.get(next(iter(col))>>2,0)+local_value(next(iter(col)),signed)) for col in cols])
        affine=np.einsum('i,ijk->jk',np.array(y[65:138],float),hop[key])+a+np.einsum('i,ijk->jk',np.array(y[:2],float)-np.array(x[2:4],float),ds)+np.diag(np.array(y[4:65],float)@ts)
        difference=max(difference,float(np.max(abs(direct+np.diag(correction)-affine))))
        ev=eigh(direct+np.diag(correction)+float(y[2])*ph+float(y[3])*q,subset_by_index=[0,0],eigvals_only=True)[0];minimum=min(minimum,(float(ev),key))
    if difference>1e-9:raise ValueError('Fresh physical matrices disagree with the search chart')
    ell=F(math.floor(minimum[0]*10**7)-1,10**7);out=dict(c,kind='hubbard_projector_extension_v17',coherent_projector={label:str(v) for label,v in zip(COHERENT,y[136:138]) if v},three_spectator_hopping={label:str(v) for label,v in zip(THREE,y[118:136]) if v},two_spectator_hopping={label:str(v) for label,v in zip(TWO,y[88:118]) if v},pair_transfer={label:str(v) for label,v in zip(PAIRS,y[84:88]) if v},spectator_hopping={label:str(v) for label,v in zip(SPECTATORS,y[70:84]) if v},hopping_telescope=str(y[65]),spin_telescope={label:str(v) for label,v in zip(LABELS,y[66:70]) if v},penalty=str(y[2]),penalized_lower=str(ell),local_window=dict(c['local_window']),joint=dict(c['joint']))
    out['joint']['penalty']=str(y[3]);out['local_window']['hopping_profile']=list(map(str,profiles(point)[1]));out['telescoping_diagonal']={str(s):str(v) for s,v in sparse.items() if v};out['signed_charge_telescope']={key:str(v) for key,v in signed.items() if v}
    lower=(ell-y[2]*theta_h-y[3]*theta_j)/5
    certificate_path=args.output/'profile_joint_r1_2_certificate.json';certificate_path.write_text(json.dumps(out,indent=2)+'\n')
    result={'accepted':False,'matrix_evaluations':calls,'optimizer_success':False if outcome is None else bool(outcome.success),'optimizer_message':'hard evaluation budget' if outcome is None else str(outcome.message),'proposed_periodic_lower':str(lower),'proposed_periodic_lower_float':float(lower),'max_fresh_affine_difference':difference,'hopping_telescope':str(y[65]),'spin_telescope':out['spin_telescope'],'spectator_hopping':out['spectator_hopping'],'pair_transfer':out['pair_transfer'],'two_spectator_hopping':out['two_spectator_hopping'],'coherent_projector':out['coherent_projector'],'three_spectator_hopping':out['three_spectator_hopping'],'signed_charge_patterns':len(out['signed_charge_telescope']),'shape_indices':indices,'seconds':time.monotonic()-started,'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name("signed_charge_numeric.py"),ROOT/"experiments/marginal_hopping_telescope.py",ROOT/"experiments/marginal_spin_telescope.py",ROOT/"experiments/marginal_spectator_hopping.py",ROOT/"experiments/marginal_pair_transfer.py",ROOT/"experiments/marginal_two_spectator_hopping.py",ROOT/"experiments/marginal_three_spectator_hopping.py",ROOT/"experiments/marginal_coherent_projector_telescope.py",args.seed,candidates_path]},'scope':'Nonaccepting numerical search in all 52 disjoint signed-charge pattern directions and nine sparse shapes. Old diagonal profiles and compact corrections are fixed, since their variation is contained in the complete signed-charge span. Two hopping profiles, both penalties, one mean-zero range-three hopping telescope four spin-dot telescopes and fourteen spectator hopping telescopes plus four pair-transfer and thirty two-spectator and eighteen three-spectator hopping telescopes plus two fixed coherent-projector telescopes are optimized. Full-spectrum soft-min annealing with signed-charge search bounds expanded to +/-8 and other original bounds retained, full spectral-Hessian coordinate scaling with explicit linear parameter bounds and finite-difference checks, at temperatures1e-6,1e-7, at most120iterations per stage and500full-spectrum evaluations total including one reserved atom export; no numerical optimality acceptance; fresh physical reconstruction; exact PSD acceptance remains separate.'}
    (args.output/'thermal_proposal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
