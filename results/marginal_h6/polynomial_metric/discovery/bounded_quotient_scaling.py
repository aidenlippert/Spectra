"""Exact coefficient-only construction and reconstruction of degree-capped quotients."""
import json,time,random
from pathlib import Path
from fractions import Fraction as F
from math import comb
from experiments.marginal_number_quotient import NumberSliceQuotient,_add
rows=[];rng=random.Random(1729)
for m,d in [(6,2),(8,3),(12,3),(16,3),(24,2),(32,2)]:
 start=time.monotonic();q=NumberSliceQuotient(list(range(m)),m//2,d);built=time.monotonic();poly={mask:F(rng.randrange(-10,11),13) for mask in rng.sample(q.monomials,min(25,len(q.monomials)))};poly={k:v for k,v in poly.items() if v};r,ideal=q.reduce(poly);reconstructed=dict(r)
 for mask,value in ideal.items():
  _add(reconstructed,{mask:F(-m//2)*value})
  for i in range(m):_add(reconstructed,{mask|(1<<i):value})
 assert reconstructed==poly;assert max(map(int.bit_count,r),default=0)<=d;assert max(map(int.bit_count,ideal),default=0)<=d-1
 rows.append(dict(q.stats,construction_seconds=built-start,reduction_and_identity_seconds=time.monotonic()-built,physical_slice_dimension_not_enumerated=comb(m,m//2),remainder_terms=len(r),witness_terms=len(ideal),exact_coefficient_identity=True))
out=Path('results/marginal_h6/polynomial_metric/bounded_quotient');out.mkdir(exist_ok=True);(out/'scaling.json').write_text(json.dumps({'scope':'One-spin bounded-degree coefficient quotients; no positivity claim. Physical-slice dimension computed symbolically, no assignments generated.','rows':rows},indent=2)+'\n');print(json.dumps(rows,indent=2))
