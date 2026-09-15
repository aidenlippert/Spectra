"""Bounded determinant-free Slater upper optimization.

The objective uses predecoded CAR words and Wick's theorem.  It never builds a
many-body vector; integer witness expansion remains confined to the separate
pair_rotation.py diagnostic.
"""
from pathlib import Path
import json, math, time
import numpy as np
from scipy.optimize import minimize
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from experiments.marginal_symbolic import decode

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'results/global_response_20260913/upper_model'

def wick(words,R):
    def pair(a,b):
        ca,ia=a; cb,ib=b
        if ca and cb or (not ca and not cb): return 0.
        if ca: return R[ib,ia]
        return (ia==ib)-R[ia,ib]
    def rec(q):
        if not q:return 1.
        if len(q)&1:return 0.
        return sum((-1)**(j-1)*pair(q[0],q[j])*rec(q[1:j]+q[j+1:]) for j in range(1,len(q)))
    return rec(words)

def prepared(data):
    # Keep words and floating coefficients immutable across objective calls.
    return [(w,float(c)) for w,c in decode(data['hamiltonian'],data['modes'],4).items()]

def objective(x,rows,modes,na,nb):
    if modes % 4 or na != nb or na != modes//4 or len(x) != modes:
        raise ValueError('This numerical proposer supports neutral closed-shell H controls only')
    ns=modes//2; R=np.zeros((modes,modes)); col=0
    for spin, n in ((0,na),(1,nb)):
        for i,t in enumerate(x[spin*ns:(spin+1)*ns]):
            if i>=n: continue
            o=2*i+spin; q=2*(i+ns//2)+spin; z=math.sqrt(1+t*t)
            R[o,o]+=1/z**2; R[q,q]+=t*t/z**2; R[o,q]+=t/z**2; R[q,o]+=t/z**2; col+=1
    return sum(c*wick(w,R) for w,c in rows)

def run(path,name,na,nb):
    data=json.loads(Path(path).read_text()); rows=prepared(data); ns=data['modes']//2
    rng=np.random.default_rng(20260913+len(name)); best=None; calls=0; start=time.monotonic()
    def f(x):
        nonlocal calls
        calls+=1; return objective(x,rows,data['modes'],na,nb)
    seeds=[np.zeros(2*ns),rng.normal(0,.05,2*ns)]
    for x in seeds:
        if time.monotonic()-start>55: break
        r=minimize(f,x,method='BFGS',options={'maxiter':40,'gtol':1e-7})
        if best is None or r.fun<best.fun: best=r
    tan=[int(round(float(v)*1e8)) for v in best.x]
    hf=f(np.zeros(2*ns)); return {'case':name,'modes':data['modes'],'particles':data['particles'],'na':na,'nb':nb,'float_energy':float(best.fun),'hf_energy':float(hf),'gain_mHa':1000*(hf-best.fun),'tangents_float':best.x.tolist(),'tangents_1e8':tan,'optimizer_success':bool(best.success),'optimizer_message':best.message,'objective_calls':calls,'seconds':time.monotonic()-start,'decoded_terms':len(rows),'construction':'paired Slater rotations from HF; no inherited amplitudes'}

if __name__=='__main__':
    cases=[('results/certificate_scaling/active_space_ladder/h6/fixture.json','h6',3,3),('results/response_consistency_20260913/fresh_h6_1p6/fixture.json','h6_1p6',3,3),('results/certificate_scaling/active_space_ladder/h8/fixture.json','h8',4,4),('results/global_response_20260913/fresh_h6_1p73/fixture.json','h6_1p73',3,3)]
    rows=[run(*c) for c in cases]; OUT.mkdir(parents=True,exist_ok=True); (OUT/'compact_trials.json').write_text(json.dumps(rows,indent=2)+'\n'); print(json.dumps(rows,indent=2))
