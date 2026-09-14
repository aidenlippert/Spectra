from itertools import combinations
from fractions import Fraction as F
from math import comb,factorial
import json
from pathlib import Path

def falling(n,k):
 result=1
 for j in range(k):result*=n-j
 return result

def average(m,p,d,a,b,c):
 return F(falling(d,c)*sum((-1)**j*comb(a-c,j)*comb(b-c,j)*factorial(j)*falling(d-c,j)*falling(p-c-j,a-c-j)*falling(p-c-j,b-c-j) for j in range(min(a-c,b-c)+1)),falling(m,a+b-c))
rows=[]
for m in [4,6,8]:
 p=m//2;subsets=[sum(1<<i for i in indices) for indices in combinations(range(m),p)];states=[(a,b,(a&b).bit_count()) for a in subsets for b in subsets];totals=[sum(d==x for a,b,d in states) for x in range(p+1)];checked=0
 for a in range(p+1):
  for b in range(p+1):
   for c in range(min(a,b)+1):
    alpha=(1<<a)-1;beta=(1<<c)-1|(((1<<(b-c))-1)<<a)
    for d in range(p+1):
     exact=F(sum(aa&alpha==alpha and bb&beta==beta and dd==d for aa,bb,dd in states),totals[d]);assert exact==average(m,p,d,a,b,c);checked+=1
 for d in range(p+1):
  value=F(1)
  for j in range(1,p+1):value*=1-F(d,j)
  assert value==(1 if d==0 else 0)
 rows.append({'sites':m,'spin_population':p,'states_enumerated':len(states),'averaging_identities_checked':checked,'projector_occupation_degree':2*p})
result={'checks':rows,'scope':'Exact finite validation of the site-permutation averaging formula and degree-m projector upper representation. General degree lower bound follows from the written polynomial-degree argument, not from extrapolation of these finite checks.'}
Path('/Users/aidenlippert/Documents/Spectra/results/marginal_h6/polynomial_metric/projector_degree_identity.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
