"""Bounded quadratic Pauli-norm diagnostic. Not an accepted dynamics checker."""
from fractions import Fraction as F
from time import perf_counter
import json
from .pauli import multiply
from .certificates import _anti,_sqrt_interval
from .v7_certificate import Generator,norm_witness
from .v7_headroom import model,ROOT,TOL


def square(op,cap=4096):
 keys=sorted(op);n=len(keys[0]) if keys else 1;out={'I'*n:sum((c*c for c in op.values()),F(0))};pairs=0
 for i,p in enumerate(keys):
  for q in keys[i+1:]:
   pairs+=1
   if _anti(p,q):continue
   phase,z=multiply(p,q)
   if phase.imag:raise AssertionError('commuting Hermitian Paulis product not real')
   out[z]=out.get(z,F(0))+2*op[p]*op[q]*int(phase.real)
   if not out[z]:del out[z]
   if len(out)>cap:raise ValueError('quadratic support budget')
 return {p:c for p,c in out.items() if c},pairs

def bound(op):
 start=perf_counter();sq,pairs=square(op)
 # Independently certify R^2 by conventional group bounds.
 options=[]
 for grouping in ('l1','firstfit','weighted'):
  w,c=norm_witness(sq,grouping);b=sum((F(x['upper']) for x in w),F(0));options.append((b,grouping,c))
 b,grouping,cost=min(options,key=lambda z:z[0]);hi=_sqrt_interval(b,16)[1]
 return hi,dict(pair_tests=pairs,square_terms=len(sq),grouping=grouping,grouping_cost=cost,elapsed=perf_counter()-start)

def run():
 source=json.loads((ROOT/'results/v7/strong_baseline.json').read_text());rows=[]
 for case in source['rows']:
  if case['status']!='certified' or case['n']>4:continue
  n,family=case['n'],case['family'];T=F(case['T']);gamma=F(case['gamma']);m=case['order']-1
  h,o=model(n,family);g=Generator(h,gamma,n,max_terms=512);current=o
  for j in range(m):current={p:c/F(j+1) for p,c in g.apply(current).items()}
  residual=g.apply(current);factor=T**(m+1)/F(m+1)
  conventional=[]
  for grouping in ('l1','firstfit','weighted'):
   w,_=norm_witness(residual,grouping);conventional.append(sum((F(x['upper']) for x in w),F(0))*factor)
  q,cost=bound(residual);q*=factor
  row=dict(n=n,family=family,gamma=str(gamma),T=str(T),order=m,terms=len(residual),best_group_bound=str(min(conventional)),quadratic_bound=str(q),
           recovers_lower_order=min(conventional)>TOL and q<=TOL,quadratic_cost=cost)
  rows.append(row)
 result=dict(rows=rows,claim='norm diagnostic only; full construction/solve/independent checking headroom NOT ESTABLISHED')
 (ROOT/'results/v7/quadratic_probe.json').write_text(json.dumps(result,indent=2));print(json.dumps({'cases':len(rows),'recovered_lower_order':sum(r['recovers_lower_order'] for r in rows),'tightened':sum(F(r['quadratic_bound'])<F(r['best_group_bound']) for r in rows),'quadratic_seconds':sum(r['quadratic_cost']['elapsed'] for r in rows)},indent=2))
if __name__=='__main__':run()
