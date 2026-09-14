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
seen=set()
for j in range(a.shape[1]):
 lo,hi=a.indptr[j:j+2];seen.add((tuple(a.indices[lo:hi]),tuple(a.data[lo:hi])))
print('distinct exact columns',len(seen),'of',a.shape[1])
