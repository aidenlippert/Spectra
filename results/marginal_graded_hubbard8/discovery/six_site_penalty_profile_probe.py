from pathlib import Path
import json,sys,math
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
BASE=ROOT/'results/marginal_graded_hubbard8'
C=json.loads((BASE/'six_site_projector/certificate.json').read_text())
vec={int(k):v for k,v in C['vector'].items()}; norm=sum(x*x for x in vec.values())
theta=1084870113/2000000000
# Local profile sums are 20, 5, 5/2 for bulk U4,t1,V1/2.

def prof(x):
 u0,u1,t0,t1,v0,v1=x
 u0,u1,t0,t1,v0,v1=[Fraction(str(z)) if not isinstance(z,Fraction) else z for z in (u0,u1,t0,t1,v0,v1)]
 return [u0,u1,10-u0-u1,10-u0-u1,u1,u0],[t0,t1,5-2*t0-2*t1,t1,t0],[v0,v1,Fraction(5,2)-2*v0-2*v1,v1,v0]

def build(x):
 op,hop,den=prof(x)
 if min(op+hop+den)<0:return None
 ss=sector_matrices(6,Fraction(10,3),1,op,hop,Fraction(1,2),den)
 out=[]
 for key,(_,K,cols) in ss.items():
  sc=np.sqrt([sum(a*a for a in c.values()) for c in cols]); A=np.asarray(K,float)/sc[:,None]/sc[None,:]
  w=np.array([sum(a*vec.get(s,0) for s,a in c.items()) for c in cols])/sc/np.sqrt(norm)
  out.append((key,A,np.outer(w,w)))
 return out
from fractions import Fraction
# precompute affine matrices by six parameter basis via building each; fit K, P fixed
base=np.array([Fraction(3,20),Fraction(83,25),Fraction(5,12),Fraction(5,4),Fraction(5,24),Fraction(5,8)],dtype=object)
points=[base]+[base+np.array([Fraction(1,100) if i==j else 0 for i in range(6)],dtype=object) for j in range(6)]
allm=[build(p) for p in points]
# P is independent; affine interpolate A
mats=[]
for j in range(len(allm[0])):
 key=allm[0][j][0]; A0=allm[0][j][1]; P=allm[0][j][2]
 As=[(allm[i+1][j][1]-A0)/0.01 for i in range(6)]
 mats.append((key,A0,As,P))
calls=0
def minimum(x,k):
 global calls;calls+=1
 if calls>120:raise ValueError('Hard evaluation cap exceeded')
 p=prof(x)
 if min(sum(p,[]))<0 or not 0<=k<=2:return -1e3,None
 best=(1e9,None)
 for key,A0,As,P in mats:
  A=A0+sum(float(x[i]-base[i])*As[i] for i in range(6))+float(k)*P
  z=float(eigh(A,subset_by_index=[0,0],eigvals_only=True)[0])
  if z<best[0]:best=(z,key)
 return best
def obj(y):
 x=y[:6];k=y[6]; ell,_=minimum(x,k);return -(ell-k*theta)/5
x0=np.array([float(z) for z in base]+[.152635])
r=minimize(obj,x0,method='Nelder-Mead',options={'maxfev':116,'xatol':1e-7,'fatol':1e-8})
baseline,baseline_sector=minimum(x0[:6],x0[6])
candidate=[Fraction(round(float(z)*10**6),10**6) for z in r.x]
ell,key=minimum(candidate[:6],candidate[6])
# Fresh exact CAR matrices, not the affine approximation, crosscheck the result.
fresh=build(candidate[:6])
direct=min((float(eigh(A+float(candidate[6])*P,subset_by_index=[0,0],eigvals_only=True)[0]),key)
           for key,A,P in fresh)
if abs(direct[0]-ell)>1e-9:raise ValueError('Fresh exact matrix differs from affine proposal')
density=(ell-float(candidate[6])*theta)/5
if density>-.6106763470511881:raise ValueError('Proposal violates existing physical upper')
out={'raw_candidate':r.x.tolist(),'rational_candidate':list(map(str,candidate)),'objective_evaluations':calls,'baseline':{'x':x0.tolist(),'minimum':baseline,'sector':baseline_sector,'periodic_density':(baseline-x0[6]*theta)/5},
     'minimum':ell,'sector':key,'periodic_density':density,'fresh_matrix_minimum':direct[0],
     'scope':'Numerical proposal only; normalized reflection blocks; fixed projector and theta; hard120 evaluation cap; fresh exact-matrix numerical crosscheck.'}
(BASE/'six_site_projector').mkdir(exist_ok=True)
(BASE/'six_site_projector/profile_probe.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out),flush=True)
