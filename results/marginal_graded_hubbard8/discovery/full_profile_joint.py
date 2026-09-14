from pathlib import Path
from fractions import Fraction as F
import json,sys,math
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import _projector_vector
from experiments.marginal_charged_projectors import charged_vectors
BASE=ROOT/'results/marginal_graded_hubbard8'; C=json.loads((BASE/'six_site_projector/refined_certificate.json').read_text())
vec={int(k):v for k,v in C['vector'].items()}; norm=sum(x*x for x in vec.values()); charged={int(k):v for k,v in json.loads((BASE/'charged_projector/source.json').read_text())['physical_state'].items()}; cvs,cn=charged_vectors(charged)
theta_h=float(F(C['projector_sum_ceiling'])/C['windows'])
def prof(x):
 a,b,p,q,d,e=x
 return [a,b,10-a-b,10-a-b,b,a],[p,q,5-2*p-2*q,q,p],[d,e,F(5,2)-2*d-2*e,e,d]
def build(x):
 op,ho,de=prof(x); cv=lambda z:F(round(float(z)*10**6),10**6); ss=sector_matrices(6,F(10,3),1,list(map(cv,op)),list(map(cv,ho)),F(1,2),list(map(cv,de)))
 out=[]
 for key,(_,K,cols) in ss.items():
  sc=np.sqrt([sum(a*a for a in c.values()) for c in cols]); A=np.asarray(K,float)/sc[:,None]/sc[None,:]
  def P(v,n):
   w=np.array([sum(a*v.get(s,0) for s,a in c.items()) for c in cols])/sc/np.sqrt(n); return np.outer(w,w)
  out.append((key,A,P(vec,norm),sum((P(v,cn) for v in cvs),np.zeros_like(A))))
 return out
base=np.array([.204105,3.039489,.478826,1.260412,.147254,.372037]); points=[base]+[base+np.eye(6)[i]*.01 for i in range(6)]; mm=[build(x) for x in points]; mats=[]
for j in range(len(mm[0])):
 k,A,P,Q=mm[0][j]; As=[(mm[i+1][j][1]-A) for i in range(6)]; mats.append((k,A,As,P,Q))
def run(ratio, baseline=False):
 # numerical source Gram ceiling from previous receipt
 rows=json.loads((BASE/'joint_projector/numeric_proposal.json').read_text())['rows']; row=next(z for z in rows if z['ratio']==str(F(ratio))); theta=float(F(row['theta_joint']))
 calls=0
 def fun(y, baseline=False):
  nonlocal calls; calls+=1
  if calls>252: raise RuntimeError('evaluation cap')
  x=y[:6]; a,b,p,q,d,e=x; k=y[6] if baseline else y[6]+y[7]; beta=0 if baseline else y[7]
  if min(prof(x)[0]+prof(x)[1]+prof(x)[2])<0 or min(y[6:])<0:return 1000
  best=1e9
  for key,A,As,P,Q in mats:
   M=A+sum((x[i]-base[i])*As[i] for i in range(6))+k*P+beta*float(ratio)*Q
   best=min(best,float(eigh(M,subset_by_index=[0,0],eigvals_only=True)[0]))
  return -(best-(y[6]*theta_h if baseline else y[6]*theta_h+y[7]*theta))/5
 x0=np.r_[base,[.152635]] if baseline else np.r_[base,[.152635,.1]]
 res=minimize(lambda y:fun(y,baseline),x0,method='Nelder-Mead',options={'maxfev':250,'xatol':1e-6,'fatol':1e-8})
 y=res.x; best=-fun(y,baseline); density=-res.fun; return {'ratio':str(ratio),'baseline':baseline,'raw':y.tolist(),'rational':[str(F(round(float(z)*10**6),10**6)) for z in y],'periodic_density':density,'evaluations':calls}
out=[run(F(1,2),True),run(F(1,2),False),run(F(1),True),run(F(1),False),run(F(2),True),run(F(2),False)]
(BASE/'joint_projector/full_profile_proposal.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
