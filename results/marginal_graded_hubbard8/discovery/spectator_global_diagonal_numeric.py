"""Nonaccepting bounded low-energy-subspace conic search; full-Fock export check."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,sys,time
import numpy as np
from scipy.linalg import eigh
import cvxpy as cp
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from joint_profile_numeric import build,profiles
from experiments.marginal_charged_projectors import charged_vectors
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value
from signed_charge_numeric import prepare
from experiments.marginal_hopping_telescope import actions,projected_matrix
from experiments.marginal_spin_telescope import LABELS,actions as spin_actions
from experiments.marginal_spectator_hopping import LABELS as SPECTATORS,actions as spectator_actions
BASE=ROOT/'results/marginal_graded_hubbard8'



def main():
    parser=argparse.ArgumentParser();parser.add_argument('seed',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();c=json.loads(args.seed.read_text())
    indices=[1,2,4,5,6,7,31,11,20]
    candidates_path=BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json';candidates=json.loads(candidates_path.read_text())['candidates']
    shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in indices]
    x,physical,mats,*_=prepare(c,shapes)
    action=[actions()]+[spin_actions({label:1}) for label in LABELS]+[spectator_actions({label:1}) for label in SPECTATORS];hop={}
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);hop[key]=np.array([projected_matrix(cols,item) for item in action],float)/scale[:,None]/scale[None,:]
    theta_h=F(c['projector_sum_ceiling'])/c['windows'];theta_j=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
    diagonal={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};taus=[]
    for shape in shapes:
        ratios={diagonal.get(s,F(0))/v for s,v in shape.items()}
        if len(ratios)!=1:raise ValueError('Seed sparse correction is outside the fixed span')
        taus.append(float(ratios.pop()))
    y0=np.array(list(map(float,x[2:4]))+[float(F(c['penalty'])),float(F(c['joint']['penalty']))]+taus+[float(F(c.get('signed_charge_telescope',{}).get(key,0))) for key in PATTERNS]+[float(F(c.get('hopping_telescope',0)))]+[float(F(c.get('spin_telescope',{}).get(label,0))) for label in LABELS]+[float(F(c.get('spectator_hopping',{}).get(label,0))) for label in SPECTATORS]+[float(F(c['penalized_lower']))])
    assert len(y0)==85
    cost=np.array([0.,0.,float(theta_h)/5,float(theta_j)/5]+[0.]*80+[-.2]);cached=None;values=None;jac=None;best=None;calls=0
    def evaluate(y):
        nonlocal cached,values,jac,best,calls
        if cached is not None and np.array_equal(y,cached):return values,jac
        calls+=1
        if calls>250:raise RuntimeError('Hard250 matrix-evaluation budget exceeded')
        cached=y.copy();vals=[];grads=[]
        for key,a,ds,ph,q,ts,*_ in mats:
            matrix=a+np.einsum('i,ijk->jk',y[:2]-np.array(x[2:4],float),ds)+y[2]*ph+y[3]*q+np.diag(y[4:65]@ts)+np.einsum('i,ijk->jk',y[65:84],hop[key])
            ev,vectors=eigh(matrix,subset_by_index=[0,0]);w=vectors[:,0]
            vals.append(float(ev[0])-y[-1]);grads.append([float(w@d@w) for d in ds]+[float(w@ph@w),float(w@q@w)]+list(ts@(w*w))+[float(w@matrix@w) for matrix in hop[key]]+[-1.])
        values=np.array(vals);jac=np.array(grads);density=(min(values)+y[-1]-y[2]*float(theta_h)-y[3]*float(theta_j))/5
        if y[0]+y[1]<=2.5+1e-12 and (best is None or density>best[0]):best=(density,y.copy())
        return values,jac
    # Solve only small compressed cones, enriching from full spectral checks.
    basis={};history=[];center=y0.copy();radius=.05
    def full_matrix(y,item):
        key,a,ds,ph,q,ts,*_=item
        return a+np.einsum('i,ijk->jk',y[:2]-np.array(x[2:4],float),ds)+y[2]*ph+y[3]*q+np.diag(y[4:65]@ts)+np.einsum('i,ijk->jk',y[65:84],hop[key])
    def enrich(y):
        spectra=[]
        for item in mats:
            matrix=full_matrix(y,item)
            ev,v=eigh(matrix,subset_by_index=[0,min(3,len(matrix)-1)])
            spectra.append((float(ev[0]),item[0],v))
        initial=not basis
        minimum=min(e for e,key,v in spectra)
        for e,key,v in spectra:
            if not initial and e>max(minimum+.01,y[-1]+1e-7):continue
            joined=np.column_stack([basis[key],v]) if key in basis else v
            u,s,_=np.linalg.svd(joined,full_matrices=False)
            n=int(sum(s>1e-8))
            if n>48:raise RuntimeError('subspace dimension cap48 reached')
            basis[key]=u[:,:n]
        if len(basis)>94 or sum(v.shape[1] for v in basis.values())>2000:raise RuntimeError('total compressed order cap2000 reached')
    evaluate(center)
    outcome=SimpleNamespace(success=False,message='bounded subspace round limit')
    for iteration in range(12):
        try:enrich(center)
        except RuntimeError as error:
            outcome.message=str(error);break
        var=cp.Variable(85)
        constraints=[var[:2]>=0,var[:2]<=2.5,cp.sum(var[:2])<=2.5,var[2:4]>=0,var[2:4]<=3,var[4:84]>=-2,var[4:84]<=2,var[-1]>=-10,var[-1]<=0]
        for item in mats:
            key,a,ds,ph,q,ts,*_=item
            diagonal_operators=list(ds)+[ph,q]+[np.diag(t) for t in ts]+list(hop[key])+[-np.eye(len(a))]
            diagonal_coefficients=np.array([np.diag(op) for op in diagonal_operators]).T
            diagonal_fixed=np.diag(a-np.einsum('i,ijk->jk',np.array(x[2:4],float),ds))
            constraints.append(diagonal_fixed+diagonal_coefficients@var>=0)
            if key not in basis:continue
            v=basis[key];n=v.shape[1]
            operators=list(ds)+[ph,q]+[np.diag(t) for t in ts]+list(hop[key])+[-np.eye(len(a))]
            coefficient=np.array([(v.T@op@v).ravel() for op in operators]).T
            fixed=v.T@(a-np.einsum('i,ijk->jk',np.array(x[2:4],float),ds))@v
            constraints.append(fixed+cp.reshape(coefficient@var,(n,n),order='C')>>0)
        problem=cp.Problem(cp.Minimize(cost@var),constraints)
        var.value=center
        problem.solve(solver='CLARABEL',time_limit=15.,max_iter=60,tol_gap_abs=1e-9,tol_feas=1e-9,verbose=False)
        entry={'round':iteration,'status':problem.status,'active_blocks':len(basis),'total_compressed_order':sum(v.shape[1] for v in basis.values()),'largest_compressed_dimension':max(v.shape[1] for v in basis.values()),'solver_seconds':problem.solver_stats.solve_time}
        if var.value is None:
            history.append(entry);outcome.message='conic solver returned no candidate';break
        center=var.value
        vals,_=evaluate(center)
        entry.update(model_periodic_upper=-float(problem.value),full_minimum_slack=float(min(vals)),best_periodic_lower=best[0])
        history.append(entry);print(json.dumps(entry),flush=True)
        if min(vals)>-1e-8:
            outcome.success=True;outcome.message='numerical full spectral feasibility at conic model candidate; not exact optimality';break
    (args.output/'subspace_basis.json').write_text(json.dumps([{'sector':key,'basis':value.tolist()} for key,value in basis.items()])+'\n')
    (args.output/'subspace_history.json').write_text(json.dumps(history,indent=2)+'\n')
    y=[F(round(float(value)*10**6),10**6) for value in best[1][:-1]]
    point=x.copy();point[2:4]=y[:2];fresh=physical(point);signed=dict(zip(PATTERNS,y[13:65]));sparse={s:sum(t*shape.get(s,0) for t,shape in zip(y[4:13],shapes)) for s in set().union(*shapes)}
    fresh_action=[actions()]+[spin_actions({label:1}) for label in LABELS]+[spectator_actions({label:1}) for label in SPECTATORS]
    minimum=(float('inf'),None);difference=0.
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        k=fresh[key][1];scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);direct=np.array(k,float)/scale[:,None]/scale[None,:]+np.einsum('i,ijk->jk',np.array(y[65:84],float),np.array([projected_matrix(cols,item) for item in fresh_action],float))/scale[:,None]/scale[None,:]
        correction=np.array([float(sparse.get(next(iter(col))&1023,0)-sparse.get(next(iter(col))>>2,0)+local_value(next(iter(col)),signed)) for col in cols])
        affine=np.einsum('i,ijk->jk',np.array(y[65:84],float),hop[key])+a+np.einsum('i,ijk->jk',np.array(y[:2],float)-np.array(x[2:4],float),ds)+np.diag(np.array(y[4:65],float)@ts)
        difference=max(difference,float(np.max(abs(direct+np.diag(correction)-affine))))
        ev=eigh(direct+np.diag(correction)+float(y[2])*ph+float(y[3])*q,subset_by_index=[0,0],eigvals_only=True)[0];minimum=min(minimum,(float(ev),key))
    if difference>1e-9:raise ValueError('Fresh physical matrices disagree with the search chart')
    ell=F(math.floor(minimum[0]*10**7)-1,10**7);out=dict(c,kind='hubbard_projector_extension_v13',spectator_hopping={label:str(v) for label,v in zip(SPECTATORS,y[70:84]) if v},hopping_telescope=str(y[65]),spin_telescope={label:str(v) for label,v in zip(LABELS,y[66:70]) if v},penalty=str(y[2]),penalized_lower=str(ell),local_window=dict(c['local_window']),joint=dict(c['joint']))
    out['joint']['penalty']=str(y[3]);out['local_window']['hopping_profile']=list(map(str,profiles(point)[1]));out['telescoping_diagonal']={str(s):str(v) for s,v in sparse.items() if v};out['signed_charge_telescope']={key:str(v) for key,v in signed.items() if v}
    lower=(ell-y[2]*theta_h-y[3]*theta_j)/5
    certificate_path=args.output/'profile_joint_r1_2_certificate.json';certificate_path.write_text(json.dumps(out,indent=2)+'\n')
    result={'accepted':False,'matrix_evaluations':calls,'optimizer_success':False if outcome is None else bool(outcome.success),'optimizer_message':'hard evaluation budget' if outcome is None else str(outcome.message),'proposed_periodic_lower':str(lower),'proposed_periodic_lower_float':float(lower),'max_fresh_affine_difference':difference,'hopping_telescope':str(y[65]),'spin_telescope':out['spin_telescope'],'spectator_hopping':out['spectator_hopping'],'signed_charge_patterns':len(out['signed_charge_telescope']),'shape_indices':indices,'seconds':time.monotonic()-started,'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name("signed_charge_numeric.py"),ROOT/"experiments/marginal_hopping_telescope.py",ROOT/"experiments/marginal_spin_telescope.py",ROOT/"experiments/marginal_spectator_hopping.py",args.seed,candidates_path]},'scope':'Nonaccepting numerical search in all 52 disjoint signed-charge pattern directions and nine sparse shapes. Old diagonal profiles and compact corrections are fixed, since their variation is contained in the complete signed-charge span. Two hopping profiles, both penalties, one mean-zero range-three hopping telescope four spin-dot telescopes and fourteen spectator hopping telescopes are optimized. At most twelve conic subspace rounds, all94blocks initialized, total compressed order2000, maximum dimension48; all4096 projected diagonal inequalities and original search bounds with no temporary trust box and 15 solver seconds/60iterations per round; fresh physical reconstruction; exact PSD acceptance remains separate.'}
    (args.output/'subspace_proposal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
