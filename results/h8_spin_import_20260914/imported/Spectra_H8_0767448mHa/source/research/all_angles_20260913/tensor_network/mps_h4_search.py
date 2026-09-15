import json, hashlib, itertools, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from scipy.optimize import minimize
from experiments.marginal_determinant_tree import DeterminantOracle
FIX=ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json'
OUT=ROOT/'results/all_angles_20260913/tensor_network'

def mps_vec(x,L,D):
    ts=[]; k=0
    for i in range(L):
        dl=1 if i==0 else D; dr=1 if i==L-1 else D
        n=2*dl*dr; ts.append(x[k:k+n].reshape(2,dl,dr)); k+=n
    v=np.zeros(1<<L)
    for s in range(1<<L):
        z=np.array([1.])
        for i in range(L): z=z@ts[i][(s>>i)&1]
        v[s]=z[0]
    return v

def main():
    fix=json.loads(FIX.read_text()); oracle=DeterminantOracle(fix); L=fix['modes']; states=[s for s in range(1<<L) if s.bit_count()==fix['particles']]; ix={s:i for i,s in enumerate(states)}
    H=np.zeros((len(states),len(states)))
    for s in states:
        for t,c in oracle.action(s).items(): H[ix[t],ix[s]]=float(c)
    sh=hashlib.sha256(FIX.read_bytes()).hexdigest(); rows=[]
    rng=np.random.default_rng(20260913)
    for D in (1,2):
        best=None
        n=sum(2*(1 if i==0 else D)*(1 if i==L-1 else D) for i in range(L))
        for restart in range(1):
            x=rng.normal(0,.2,n)
            def fun(y):
                w=mps_vec(y,L,D); w=w[states]; q=np.linalg.norm(w)
                if q<1e-14:return 1e6
                return float(w@(H@w)/(q*q))
            r=minimize(fun,x,method='BFGS',options={'maxiter':40,'gtol':1e-6})
            if best is None or r.fun<best.fun: best=r
        w=mps_vec(best.x,L,D); q=np.linalg.norm(w[states]); wf=w[states]/q
        den=10**10; amps=np.rint(wf*den).astype(np.int64); nz=np.flatnonzero(amps); rat=np.zeros_like(wf); rat[nz]=amps[nz]/den; rat/=np.linalg.norm(rat)
        exact=float(rat@(H@rat)); witness={'states':[states[i] for i in nz],'amplitudes':[int(a) for a in amps[nz]]}; oracle_exact=oracle.upper(witness)
        rows.append({'bond_dimension':D,'params':n,'float_energy':float(best.fun),'rounded_rational_energy':exact,'exact_integer_witness_energy':str(oracle_exact),'support':int(len(nz)),'optimizer_success':bool(best.success),'optimizer_message':best.message})
        if D==2: (OUT/'mps_h4_D2_witness.json').write_text(json.dumps(witness,indent=2)+'\n')
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'mps_h4_report.json').write_text(json.dumps({'fixture':str(FIX),'fixture_sha256':sh,'method':'random-start variational open-boundary MPS; no FCI seed; fixed-N projection by determinant filtering','rows':rows},indent=2)+'\n')
    print(json.dumps(rows,indent=2))
if __name__=='__main__': main()
