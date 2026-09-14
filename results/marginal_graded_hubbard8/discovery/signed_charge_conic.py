"""Bounded conic cross-check; physical exact acceptance remains separate."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys,time,math,hashlib
import numpy as np
import cvxpy as cp
from scipy.linalg import eigh
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from signed_charge_numeric import prepare,profiles,BASE
from experiments.marginal_signed_charge_telescope import PATTERNS,local_value


def main():
    parser=argparse.ArgumentParser();parser.add_argument('seed',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True);started=time.monotonic();c=json.loads(args.seed.read_text())
    candidates=json.loads((BASE/'joint_projector/signed_density/symmetry_diagonal_candidates.json').read_text())['candidates'];shapes=[{int(s):v for s,v in candidates[i]['diagonal'].items()} for i in [1,2,4,5,6,7,31,11,20]]
    x,physical,mats,*_=prepare(c,shapes);th=F(c['projector_sum_ceiling'])/c['windows'];tj=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows'];diagonal={int(s):F(v) for s,v in c['telescoping_diagonal'].items()}
    taus=[]
    for shape in shapes:
        ratios={diagonal.get(s,F(0))/v for s,v in shape.items()}
        if len(ratios)!=1:raise ValueError('Sparse span mismatch')
        taus.append(float(ratios.pop()))
    seed=np.array(list(map(float,x[2:4]))+[float(F(c['penalty'])),float(F(c['joint']['penalty']))]+taus+[float(F(c['signed_charge_telescope'].get(k,0))) for k in PATTERNS]+[float(F(c['penalized_lower']))])
    def evaluate(y):
        values={}
        for key,a,ds,ph,q,ts,*_ in mats:
            m=a+(y[0]-float(x[2]))*ds[0]+(y[1]-float(x[3]))*ds[1]+y[2]*ph+y[3]*q+np.diag(y[4:-1]@ts)
            values[key]=float(eigh(m,subset_by_index=[0,0],eigvals_only=True)[0])
        return values,(min(values.values())-y[2]*float(th)-y[3]*float(tj))/5
    initial,score=evaluate(seed);best=(score,seed.copy());active={key for key,v in initial.items() if v<=float(F(c['penalized_lower']))+1e-4};rounds=[]
    for step in range(2):
        if len(active)>16:raise ValueError('Bounded conic active-set size exceeded')
        y=cp.Variable(66);constraints=[y[:2]>=0,y[:2]<=2.5,cp.sum(y[:2])<=2.5,y[2:4]>=0,y[2:4]<=3,y[4:-1]>=-2,y[4:-1]<=2,y[-1]>=-10,y[-1]<=0]
        for key,a,ds,ph,q,ts,*_ in mats:
            if key in active:constraints.append(a+(y[0]-float(x[2]))*ds[0]+(y[1]-float(x[3]))*ds[1]+y[2]*ph+y[3]*q+cp.diag(ts.T@y[4:-1])-y[-1]*np.eye(len(a))>>0)
        problem=cp.Problem(cp.Minimize((y[2]*float(th)+y[3]*float(tj)-y[-1])/5),constraints);y.value=best[1]
        problem.solve(solver='CLARABEL',time_limit=40.,max_iter=60,tol_gap_abs=1e-8,tol_feas=1e-8,verbose=False)
        receipt={'round':step,'status':problem.status,'active_blocks':len(active),'solver_seconds':problem.solver_stats.solve_time}
        if y.value is None:rounds.append(receipt);break
        values,bound=evaluate(y.value);receipt.update(fresh_periodic_lower=bound,minimum_slack=min(values.values())-float(y.value[-1]));rounds.append(receipt);print(json.dumps(receipt),flush=True)
        if bound>best[0] and min(y.value[:4])>=-1e-10 and sum(y.value[:2])<=2.5+1e-10:best=(bound,y.value.copy())
        extra={key for key,v in values.items() if key not in active and v<float(y.value[-1])-1e-7}
        if not extra:break
        active|=extra
    y=[F(round(float(v)*10**6),10**6) for v in best[1][:-1]];point=x.copy();point[2:4]=y[:2];fresh=physical(point);signed=dict(zip(PATTERNS,y[13:]));sparse={s:sum(t*shape.get(s,0) for t,shape in zip(y[4:13],shapes)) for s in set().union(*shapes)};minimum=float('inf');difference=0.
    for key,a,ds,ph,q,ts,cols,*_ in mats:
        k=fresh[key][1];scale=np.sqrt([sum(v*v for v in col.values()) for col in cols]);direct=np.array(k,float)/scale[:,None]/scale[None,:];correction=np.array([float(sparse.get(next(iter(col))&1023,0)-sparse.get(next(iter(col))>>2,0)+local_value(next(iter(col)),signed)) for col in cols]);affine=a+(float(y[0])-float(x[2]))*ds[0]+(float(y[1])-float(x[3]))*ds[1]+np.diag(np.array(y[4:],float)@ts);difference=max(difference,float(np.max(abs(direct+np.diag(correction)-affine))));minimum=min(minimum,float(eigh(direct+np.diag(correction)+float(y[2])*ph+float(y[3])*q,subset_by_index=[0,0],eigvals_only=True)[0]))
    if difference>1e-9:raise ValueError('Fresh physical reconstruction mismatch')
    ell=F(math.floor(minimum*10**7)-1,10**7);out=dict(c,penalty=str(y[2]),penalized_lower=str(ell),local_window=dict(c['local_window']),joint=dict(c['joint']));out['joint']['penalty']=str(y[3]);out['local_window']['hopping_profile']=list(map(str,profiles(point)[1]));out['telescoping_diagonal']={str(s):str(v) for s,v in sparse.items() if v};out['signed_charge_telescope']={k:str(v) for k,v in signed.items() if v}
    (args.output/'profile_joint_r1_2_certificate.json').write_text(json.dumps(out,indent=2)+'\n');result={'accepted':False,'solver':'CLARABEL','rounds':rounds,'proposed_periodic_lower':str((ell-y[2]*th-y[3]*tj)/5),'proposed_periodic_lower_float':float((ell-y[2]*th-y[3]*tj)/5),'max_fresh_affine_difference':difference,'seconds':time.monotonic()-started,'source_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),args.seed]},'scope':'Bounded conic numerical cross-check, at most two rounds of40solverseconds/60iterations and16active local blocks. Final fresh eigenvalue check covers all94blocks. This output is a proposal, not exact PSD or optimality acceptance.'};(args.output/'conic_proposal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
