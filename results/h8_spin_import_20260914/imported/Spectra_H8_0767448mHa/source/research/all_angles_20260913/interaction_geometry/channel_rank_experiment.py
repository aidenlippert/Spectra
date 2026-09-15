"""Density/Coulomb channel and cut-rank measurements on H4/H6/H8 fixtures.

This is a proposal experiment: rank tails are numerical.  The rationalized
reconstruction residual is an outward-safe coefficient L1 bound for the
two-body tensor (not a claim about SOS Gram rank).
"""
import json, glob
from pathlib import Path
from fractions import Fraction
import numpy as np

ROOT=Path('/Users/aidenlippert/Documents/Spectra')
FIX=[str(ROOT/f'results/certificate_scaling/active_space_ladder/h{i}/fixture.json') for i in (4,6,8)]
FIX += [str(ROOT/f'results/certificate_scaling/active_space_ladder_boys/h{i}/fixture.json') for i in (4,6)]

def frac(x): return Fraction(x)
def tensor(d):
    m=d['modes']; T=np.zeros((m*m,m*m)); exact={}
    for row in d['hamiltonian']:
        w=row['word']
        if len(w)!=4 or any(x[0]!=1 for x in w[:2]) or any(x[0]!=0 for x in w[2:]): continue
        p,q=w[0][1],w[1][1]; r,s=w[2][1],w[3][1]
        v=float(frac(row['coefficient']))
        T[p*m+q,r*m+s]+=v
        exact[(p,q,r,s)]=exact.get((p,q,r,s),Fraction(0))+frac(row['coefficient'])
    return T,exact
def report_one(path):
    d=json.load(open(path)); T,E=tensor(d); m=d['modes']; u,s,v=np.linalg.svd(T,full_matrices=False)
    ranks=[]
    for k in [1,2,4,8, min(16,m*m),m*m]:
        k=min(k,len(s)); tail=float(np.linalg.norm(s[k:]))
        # rationalized rank-k factors, then exact coefficient L1 residual
        A=u[:,:k]*np.sqrt(s[:k]); B=(np.sqrt(s[:k])[:,None]*v[:k,:])
        den=10**8; Ar=np.rint(A*den).astype(np.int64); Br=np.rint(B*den).astype(np.int64)
        l1=Fraction(0)
        for p in range(m):
          for q in range(m):
           for r in range(m):
            for z in range(m):
             pred=sum(Fraction(int(Ar[p*m+q,a])*int(Br[a,r*m+z]),den*den) for a in range(k))
             l1 += abs(E.get((p,q,r,z),Fraction(0))-pred)
        ranks.append({'rank':k,'spectral_tail':tail,'rational_denominator':den,'tensor_l1_bound':float(l1),'tensor_l1_bound_exact':str(l1)})
    # Adverse control: a random rank-one channel with matched Frobenius norm.
    rng=np.random.default_rng(991+m); aa=rng.normal(size=m*m); bb=rng.normal(size=m*m)
    aa*=np.sqrt(np.linalg.norm(T,'fro')/(np.linalg.norm(aa)*np.linalg.norm(bb))); bb*=np.sqrt(np.linalg.norm(T,'fro')/(np.linalg.norm(aa)*np.linalg.norm(bb)))
    ar=np.rint(aa*10**8).astype(np.int64); br=np.rint(bb*10**8).astype(np.int64)
    adverse=Fraction(0)
    for p in range(m):
      for q in range(m):
       for r in range(m):
        for z in range(m): adverse += abs(E.get((p,q,r,z),Fraction(0))-Fraction(int(ar[p*m+q])*int(br[r*m+z]),10**16))
    return {'fixture':path.replace(str(ROOT)+'/',''),'modes':m,'quartic_nonzero':len(E),'frobenius':float(np.linalg.norm(T)),'singular_values':[float(x) for x in s[:min(12,len(s))]],'ranks':ranks,'adverse_random_rank1_tensor_l1':str(adverse),'basis':d.get('orbital_basis'),'localized_basis':d.get('localized_basis')}
out=ROOT/'results/all_angles_20260913/interaction_geometry'; out.mkdir(parents=True,exist_ok=True)
res=[report_one(p) for p in FIX]; (out/'channel_rank_metrics.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
