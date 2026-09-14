import json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy.sparse import csc_matrix
from scipy.linalg import qr
import cvxpy as cp
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_spin_constructor import spin_states
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');out=root/'full_sos';out.mkdir(exist_ok=True)
c=json.loads((root/'candidate.json').read_text());o=JointPolynomial(c);compiled=json.loads((root/'compiled.json').read_text());states=spin_states(o.o.oracle);sa=np.array(states,dtype=np.int64);d=np.array([sum(((s>>(2*i))&3)==3 for i in range(o.sites))-1 for s in states]);labels=[];ri=[];ci=[];values=[];start=time.monotonic()
for degree in range(7):
 for bits in combinations(range(o.modes),degree):
  mask=sum(1<<i for i in bits)
  for assignment in range(1<<degree):
   occupied=sum(1<<i for j,i in enumerate(bits) if assignment&(1<<j));active=np.flatnonzero((sa&mask)==occupied)
   if len(active)>1:
    j=len(labels);labels.append(['positive',mask,occupied]);ri.extend(active);ci.extend([j]*len(active));values.extend([1]*len(active))
   if degree<=4:
    nonzero=active[d[active]!=0]
    if len(nonzero):
     j=len(labels);labels.append(['charge',mask,occupied]);ri.extend(nonzero);ci.extend([j]*len(nonzero));values.extend(d[nonzero].tolist())
a=csc_matrix((np.array(values,dtype=float),(ri,ci)),shape=(len(states),len(labels)));exact={m:int(v) for m,v in compiled['numerator']};scale=compiled['receipt']['numerator_scale'];p=np.array([float(F(sum(v for m,v in exact.items() if m&s==m),scale)) for s in states]);print(json.dumps({'phase':'matrix','shape':a.shape,'nnz':a.nnz,'seconds':time.monotonic()-start}),flush=True)
for degree in [2,3]:
 basis=[sum(1<<i for i in inds) for k in range(degree+1) for inds in combinations(range(o.modes),k)];full=np.array([[(s&m)==m for m in basis] for s in states],dtype=float);_,r,pivot=qr(full,mode='economic',pivoting=True);rank=int(np.sum(abs(np.diag(r))>1e-8));selected=sorted(pivot[:rank]);f=full[:,selected];basis=[basis[i] for i in selected]
 # Dual signed occupation functional; all candidate atoms and all polynomial squares are constrained.
 z=cp.Variable(len(states));nonnegative=a.T@z>=0;normal=cp.sum(z)==1;psd=f.T@cp.diag(z)@f>>0;prob=cp.Problem(cp.Minimize(p@z),[normal,nonnegative,psd]);start=time.monotonic();print(json.dumps({'phase':'solve','degree':degree,'rank':rank}),flush=True)
 try:prob.solve(solver='CLARABEL',max_iter=100,time_limit=90,tol_gap_abs=1e-7,tol_feas=1e-8,tol_gap_rel=1e-7)
 except Exception as exc:
  (out/f'degree{degree}_error.json').write_text(json.dumps({'error':str(exc),'seconds':time.monotonic()-start},indent=2)+'\n');print(str(exc),flush=True);continue
 receipt={'degree':degree,'rank':rank,'status':prob.status,'value':prob.value,'seconds':time.monotonic()-start,'scope':'Numerical full SOS plus nonsingleton indicator/localizer dual on explicit400-sector coordinates. QR rank selection is numerical discovery only.'}
 if z.value is not None:
  moment=f.T@(z.value[:,None]*f);receipt.update(minimum_atom=float(np.min(a.T@z.value)),normalization=float(sum(z.value)),minimum_moment_eigenvalue=float(np.linalg.eigvalsh(moment)[0]))
  proposal={'states':states,'basis':basis,'dual_weights':z.value.tolist(),'bound':float(-normal.dual_value),'atom_labels':labels,'atom_weights':nonnegative.dual_value.tolist(),'gram':psd.dual_value.tolist()}
  (out/f'degree{degree}_proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
 (out/f'degree{degree}_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
 if prob.value is not None and prob.value>1e-5:break
