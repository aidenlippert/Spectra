"""Find a local density-matrix witness capping the fixed-projector relaxation.

Numerical LP is proposal-only. Any resulting exact rational mixture must still
be replayed independently before a family-limit claim is accepted.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_charged_projectors import charged_vectors
BASE=ROOT/'results/marginal_graded_hubbard8';OUT=BASE/'joint_projector/signed_density'


def profiles(x):
 a,b,p,q,d,e=x
 return [a,b,10-a-b,10-a-b,b,a],[p,q,5-2*p-2*q,q,p],[d,e,F(5,2)-2*d-2*e,e,d]


def build(x):
 u,t,v=profiles(x)
 return sector_matrices(6,F(10,3),1,u,t,F(1,2),v)


def solve_exact(A,b):
 n=len(b);M=[list(map(F,row))+[F(v)] for row,v in zip(A,b)]
 for j in range(n):
  pivot=next((i for i in range(j,n) if M[i][j]),None)
  if pivot is None:raise ValueError('Singular rational basis')
  M[j],M[pivot]=M[pivot],M[j];q=M[j][j];M[j]=[v/q for v in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[a-q*b for a,b in zip(M[i],M[j])]
 return [row[-1] for row in M]


def main():
 c=json.loads((OUT/'profile_joint_r1_2_certificate.json').read_text());local=c['local_window']
 half={int(k):a for k,a in c['vector'].items()};hn=sum(a*a for a in half.values());cv,cn=charged_vectors(c['joint']['vector'])
 x=list(map(F,[local['onsite_profile'][0],local['onsite_profile'][1],local['hopping_profile'][0],local['hopping_profile'][1],local['density_profile'][0],local['density_profile'][1]]))
 original=build(x);pert=[]
 for j in range(6):
  y=x.copy();y[j]+=F(1,100);pert.append(build(y))
 alpha=F(c['penalty']);beta=F(c['joint']['penalty']);ratio=F(c['joint']['ratio'])
 states=[];rows=[];energies=[]
 for key,(_,K,cols) in original.items():
  sc=np.sqrt([sum(a*a for a in col.values()) for col in cols]);A=np.array(K,float)/sc[:,None]/sc[None,:]
  exact_D=[[[100*(v-K[i][k]) for k,v in enumerate(row)] for i,row in enumerate(ss[key][1])] for ss in pert]
  Ds=np.array([np.array(D,float)/sc[:,None]/sc[None,:] for D in exact_D])
  def proj(v,n):
   w=np.array([sum(a*v.get(s,0) for s,a in col.items()) for col in cols],float)/sc/np.sqrt(float(n));return np.outer(w,w)
  P=proj(half,hn);Q=sum((proj(v,cn) for v in cv),np.zeros_like(A));M=A+float(alpha+beta)*P+float(beta*ratio)*Q
  lowest=float(eigh(M,eigvals_only=True,subset_by_index=[0,0])[0])
  if lowest>float(F(c['penalized_lower']))+1e-4:continue
  # Symmetry-equivalent sectors duplicate rows; retain their physical vectors
  # initially, then remove duplicate numerical expectation columns in the LP.
  for delta in [np.zeros(6)]+[sign*np.eye(6)[j]*.002 for j in range(6) for sign in (-1,1)]:
   _,vs=eigh(M+np.einsum('i,ijk->jk',delta,Ds),subset_by_index=[0,0]);z=vs[:,0]/sc;z=z/max(abs(z))
   coeff=[round(float(a)*10**8) for a in z];vector={}
   for a,col in zip(coeff,cols):
    for s,b in col.items():
     if a:vector[str(s)]=a*b
   norm=sum(a*a for a in vector.values())
   def exp(B):return sum((F(a*b)*B[i][j] for i,a in enumerate(coeff) if a for j,b in enumerate(coeff) if b),F(0))/norm
   deriv=[exp(D) for D in exact_D];energy=exp(K)-sum(a*b for a,b in zip(x,deriv))
   def fidelity(v,n):return F(sum(a*v.get(int(s),0) for s,a in vector.items())**2,norm*n)
   ph=fidelity(half,hn);q=ph+ratio*sum((fidelity(v,cn) for v in cv),F(0))
   col=[F(1)]+deriv+[ph,q]
   if any(max(abs(float(a-b)) for a,b in zip(col,old))<1e-12 for old in rows):continue
   states.append(vector);rows.append(col);energies.append(energy)
 print('state_count',len(states),flush=True)
 theta_h=F(c['projector_sum_ceiling'])/c['windows'];theta_j=F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
 arr=np.array(rows,float).T
 res=linprog(np.array(energies,float),A_eq=arr[:7],b_eq=[1]+[0]*6,A_ub=arr[7:],b_ub=[float(theta_h),float(theta_j)],bounds=(0,None),method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
 if not res.success:
  (OUT/'family_limit_diagnostic.json').write_text(json.dumps({'accepted':False,'status':res.message,'states':len(states)},indent=2)+'\n');print(res.message);return
 # Add inequality slack columns; reconstruct the selected LP basis exactly.
 columns=rows+[[F(0)]*7+[F(1),F(0)],[F(0)]*8+[F(1)]]
 values=list(res.x)+list(np.array([float(theta_h),float(theta_j)])-arr[7:]@res.x)
 support=[i for i,v in enumerate(values) if v>1e-10]
 diagnostic={'accepted':False,'states':len(states),'numerical_upper_density':float(res.fun)/5,'basis_size':len(support),'support':support}
 (OUT/'family_limit_diagnostic.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
 if len(support)!=9:print(json.dumps(diagnostic));return
 weights=solve_exact([[columns[j][i] for j in support] for i in range(9)],[1]+[0]*6+[theta_h,theta_j])
 if min(weights)<0:raise ValueError('Exact basis contains negative weight/slack')
 mixture=[{'weight':str(w),'vector':states[j]} for j,w in zip(support,weights) if j<len(states)]
 upper=sum(w*energies[j] for j,w in zip(support,weights) if j<len(states))/5
 proposal={'kind':'joint_projector_family_limit_v1','half_vector':c['vector'],'charged_vector':c['joint']['vector'],'ratio':str(ratio),'theta_half':str(theta_h),'theta_joint':str(theta_j),'mixture':mixture,'proposed_periodic_family_upper':str(upper),'scope':'Exact rational basis proposed from a numerical LP. Independent physical-expectation replay remains required.'}
 (OUT/'family_limit_certificate.json').write_text(json.dumps(proposal,indent=2)+'\n')
 print(json.dumps(dict(diagnostic,rational_upper_float=float(upper),weights=len(mixture))),flush=True)


if __name__=='__main__':main()
