import json,time
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy.sparse import csc_matrix
from scipy.optimize import linprog
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_spin_constructor import spin_states
root=Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric');out=root/'localized_collective_square_quotient';out.mkdir(exist_ok=True)
c=json.loads((root/'candidate.json').read_text());o=JointPolynomial(c);compiled=json.loads((root/'compiled.json').read_text());states=spin_states(o.o.oracle);sa=np.array(states,dtype=np.int64);d=np.array([sum(((s>>(2*i))&3)==3 for i in range(o.sites))-1 for s in states]);labels=[['b']];ri=list(range(len(states)));ci=[0]*len(states);values=[1]*len(states);start=time.monotonic()
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
for sites in combinations(range(o.sites),3):
 for family,center in [('doublon',2),('number',5)]:
  evaluated=np.array([sum((((s>>(2*i))&3)==3) if family=='doublon' else ((s>>(2*i))&3).bit_count() for i in sites)-center for s in states]);active=np.flatnonzero(evaluated)
  for local_degree in range(3):
   for local_bits in combinations(range(o.modes),local_degree):
    required=sum(1<<i for i in local_bits)
    for local_assignment in range(1<<local_degree):
     occupied=sum(1<<i for j,i in enumerate(local_bits) if local_assignment&(1<<j));chosen=active[(sa[active]&required)==occupied]
     if len(chosen):
      j=len(labels);labels.append(['localized_collective_square',family,list(sites),center,required,occupied]);ri.extend(chosen);ci.extend([j]*len(chosen));values.extend((evaluated[chosen]**2).tolist())
a=csc_matrix((np.array(values,dtype=float),(ri,ci)),shape=(len(states),len(labels)));exact={m:int(v) for m,v in compiled['numerator']};scale=compiled['receipt']['numerator_scale'];p=[sum(v for m,v in exact.items() if m&s==m) for s in states];objective=np.zeros(len(labels));objective[0]=-1
print(json.dumps({'phase':'matrix','rows':len(states),'columns':len(labels),'nnz':a.nnz,'seconds':time.monotonic()-start}),flush=True);start=time.monotonic();r=linprog(objective,A_eq=a,b_eq=np.array([float(F(v,scale)) for v in p]),bounds=[(None,None)]+[(0,None)]*(len(labels)-1),method='highs-ipm',options={'time_limit':60})
receipt={'success':r.success,'status':r.message,'seconds':time.monotonic()-start,'rows':len(states),'columns':len(labels),'nnz':a.nnz,'scope':'Finite explicit400-state numerical quotient; b is FREE and objective maximizes b. Nonsingleton indicators<=6, charge localizers<=4 plus collective doublon and number squares on every triple of spatial sites times all indicators of degree<=2. Numerical sector-coordinate discovery.'}
if r.success:
 receipt['optimal_b_float']=float(r.x[0]);z=[F(float(-v)).limit_denominator(10**6) for v in r.eqlin.marginals];mom={}
 for s,v in zip(states,z):
  sub=s
  while True:
   mom[sub]=mom.get(sub,F(0))+v
   if sub==0:break
   sub=(sub-1)&s
 minimum=F(0);bad=0
 for j in range(1,a.shape[1]):
  lo,hi=a.indptr[j:j+2];v=sum((z[int(i)]*int(w) for i,w in zip(a.indices[lo:hi],a.data[lo:hi])),F(0));minimum=min(minimum,v);bad+=v<0
 value=sum((zz*pp for zz,pp in zip(z,p)),F(0))/scale;receipt.update(normalization=str(sum(z,F(0))),minimum_atom=str(minimum),violations=bad,exact_value=str(value),exact_value_float=float(value),max_abs_moment=str(max(map(abs,mom.values()))))
 (out/'dual_proposal.json').write_text(json.dumps({'states':states,'weights':[str(x) for x in z]},indent=2)+'\n')
 (out/'primal_proposal.json').write_text(json.dumps({'labels':[l for l,x in zip(labels,r.x) if abs(x)>1e-12],'values':[float(x) for x in r.x if abs(x)>1e-12]},indent=2)+'\n')
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
