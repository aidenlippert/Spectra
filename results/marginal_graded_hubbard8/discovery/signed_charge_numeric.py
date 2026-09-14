"""Bounded numerical search in disjoint signed-charge pattern coordinates."""
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
BASE=ROOT/'results/marginal_graded_hubbard8'


def prepare(c,shapes):
    local=c['local_window'];x=[F(local[name][i]) for name,i in [('onsite_profile',0),('onsite_profile',1),('hopping_profile',0),('hopping_profile',1),('density_profile',0),('density_profile',1)]]
    if any(F(c['target'][key])!=value for key,value in [('U',F(4)),('t',F(1)),('V',F(1,2))]):raise ValueError('Fixed U4,t1,V1/2 required')
    def physical(point):
        return build(point,W=F(c['target']['W']),range_profile=list(map(F,local['range_two_density_profile'])),quadratic_terms={k:F(v) for k,v in c['quadratic_charge_telescope'].items()},square_terms={k:F(v) for k,v in c['charge_square_pair_telescope'].items()},indicator_terms={k:F(v) for k,v in c['higher_charge_indicator_telescope'].items()})
    original=physical(x);perturbed=[]
    for j in (2,3):
        point=x.copy();point[j]+=F(1,100);perturbed.append(physical(point))
    half={int(s):a for s,a in c['vector'].items()};hn=sum(a*a for a in half.values());charged,cn=charged_vectors(c['joint']['vector']);ratio=F(c['joint']['ratio'])
    mats=[]
    for key,(_,k,cols) in original.items():
        scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);denom=scale[:,None]*scale[None,:]
        ds=[[[100*(v-k[i][j]) for j,v in enumerate(row)] for i,row in enumerate(p[key][1])] for p in perturbed]
        ts=np.array([[shape.get(next(iter(col))&1023,0)-shape.get(next(iter(col))>>2,0) for col in cols] for shape in shapes]+[[int(local_value(next(iter(col)),{label:F(1)})) for col in cols] for label in PATTERNS])
        def projector(vector,norm):
            w=np.array([sum(a*vector.get(s,0) for s,a in col.items()) for col in cols],float)/scale/np.sqrt(float(norm));return np.outer(w,w)
        ph=projector(half,hn);q=ph+float(ratio)*sum((projector(v,cn) for v in charged),np.zeros_like(denom))
        mats.append((key,np.array(k,float)/denom,np.array(ds,float)/denom,ph,q,ts,cols,k,ds))
    return x,physical,mats,half,hn,charged,cn,ratio


def main():
    parser=argparse.ArgumentParser();parser.add_argument('seed',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    started=time.monotonic();c=json.loads(args.seed.read_text())
    indices=[1,2,4,5,6,7,31,11,20]
    candidates_path=BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json';candidates=json.loads(candidates_path.read_text())['candidates']
    shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in indices]
    x,physical,mats,*_=prepare(c,shapes)
    theta_h=F(c['projector_sum_ceiling'])/c['windows'];theta_j=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
    diagonal={int(s):F(v) for s,v in c['telescoping_diagonal'].items()};taus=[]
    for shape in shapes:
        ratios={diagonal.get(s,F(0))/v for s,v in shape.items()}
        if len(ratios)!=1:raise ValueError('Seed sparse correction is outside the fixed span')
        taus.append(float(ratios.pop()))
    y0=np.array(list(map(float,x[2:4]))+[float(F(c['penalty'])),float(F(c['joint']['penalty']))]+taus+[float(F(c.get('signed_charge_telescope',{}).get(key,0))) for key in PATTERNS]+[float(F(c['penalized_lower']))])
    assert len(y0)==66
    cost=np.array([0.,0.,float(theta_h)/5,float(theta_j)/5]+[0.]*61+[-.2]);cached=None;values=None;jac=None;best=None;calls=0
    def evaluate(y):
        nonlocal cached,values,jac,best,calls
        if cached is not None and np.array_equal(y,cached):return values,jac
        calls+=1
        if calls>250:raise RuntimeError('Hard250 matrix-evaluation budget exceeded')
        cached=y.copy();vals=[];grads=[]
        for key,a,ds,ph,q,ts,*_ in mats:
            matrix=a+np.einsum('i,ijk->jk',y[:2]-np.array(x[2:4],float),ds)+y[2]*ph+y[3]*q+np.diag(y[4:-1]@ts)
            ev,vectors=eigh(matrix,subset_by_index=[0,0]);w=vectors[:,0]
            vals.append(float(ev[0])-y[-1]);grads.append([float(w@d@w) for d in ds]+[float(w@ph@w),float(w@q@w)]+list(ts@(w*w))+[-1.])
        values=np.array(vals);jac=np.array(grads);density=(min(values)+y[-1]-y[2]*float(theta_h)-y[3]*float(theta_j))/5
        if y[0]+y[1]<=2.5+1e-12 and (best is None or density>best[0]):best=(density,y.copy())
        return values,jac
    outcome=None
    try:
        outcome=minimize(lambda y:float(cost@y),y0,jac=lambda y:cost,method='SLSQP',bounds=[(0,2.5)]*2+[(0,3)]*2+[(-2,2)]*61+[(-10,0)],constraints=[{'type':'ineq','fun':lambda y:evaluate(y)[0],'jac':lambda y:evaluate(y)[1]},{'type':'ineq','fun':lambda y:2.5-y[0]-y[1],'jac':lambda y:np.r_[-1.,-1.,np.zeros(64)]}],options={'maxiter':240,'ftol':1e-10})
    except RuntimeError as exc:
        if 'budget' not in str(exc):raise
    y=[F(round(float(value)*10**6),10**6) for value in best[1][:-1]]
    point=x.copy();point[2:4]=y[:2];fresh=physical(point);signed=dict(zip(PATTERNS,y[13:]));sparse={s:sum(t*shape.get(s,0) for t,shape in zip(y[4:13],shapes)) for s in set().union(*shapes)}
    minimum=(float('inf'),None);difference=0.
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        k=fresh[key][1];scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);direct=np.array(k,float)/scale[:,None]/scale[None,:]
        correction=np.array([float(sparse.get(next(iter(col))&1023,0)-sparse.get(next(iter(col))>>2,0)+local_value(next(iter(col)),signed)) for col in cols])
        affine=a+np.einsum('i,ijk->jk',np.array(y[:2],float)-np.array(x[2:4],float),ds)+np.diag(np.array(y[4:],float)@ts)
        difference=max(difference,float(np.max(abs(direct+np.diag(correction)-affine))))
        ev=eigh(direct+np.diag(correction)+float(y[2])*ph+float(y[3])*q,subset_by_index=[0,0],eigvals_only=True)[0];minimum=min(minimum,(float(ev),key))
    if difference>1e-9:raise ValueError('Fresh physical matrices disagree with the search chart')
    ell=F(math.floor(minimum[0]*10**7)-1,10**7);out=dict(c,kind='hubbard_projector_extension_v10',penalty=str(y[2]),penalized_lower=str(ell),local_window=dict(c['local_window']),joint=dict(c['joint']))
    out['joint']['penalty']=str(y[3]);out['local_window']['hopping_profile']=list(map(str,profiles(point)[1]));out['telescoping_diagonal']={str(s):str(v) for s,v in sparse.items() if v};out['signed_charge_telescope']={key:str(v) for key,v in signed.items() if v}
    lower=(ell-y[2]*theta_h-y[3]*theta_j)/5
    certificate_path=args.output/'profile_joint_r1_2_certificate.json';certificate_path.write_text(json.dumps(out,indent=2)+'\n')
    result={'accepted':False,'matrix_evaluations':calls,'optimizer_success':False if outcome is None else bool(outcome.success),'optimizer_message':'hard evaluation budget' if outcome is None else str(outcome.message),'proposed_periodic_lower':str(lower),'proposed_periodic_lower_float':float(lower),'max_fresh_affine_difference':difference,'signed_charge_patterns':len(out['signed_charge_telescope']),'shape_indices':indices,'seconds':time.monotonic()-started,'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),args.seed,candidates_path]},'scope':'Nonaccepting numerical search in all 52 disjoint signed-charge pattern directions and nine sparse shapes. Old diagonal profiles and compact corrections are fixed, since their variation is contained in the complete signed-charge span. Two hopping profiles and both penalties are optimized. Hard250 matrix evaluations; fresh physical reconstruction; exact PSD acceptance remains separate.'}
    (args.output/'signed_charge_proposal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
