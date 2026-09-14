"""Exact action-free polynomial degree diagnostic for local Hubbard chains."""
import json,time
from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch
from experiments.marginal_symbolic import mono,product,scale,add,encode
from experiments.marginal_polynomial_metric import JointPolynomial
rows=[]
for metric in ('D','1'):
 for sites in ([4,6,8,10] if metric=='D' else [4,6,8]):
  modes=2*sites;n=[mono(((1,i),(0,i))) for i in range(modes)]
  h=add(*(scale(product(n[2*i],n[2*i+1]),4) for i in range(sites)),*(mono(((1,2*i+s),(0,2*(i+1)+s)),F(-1,3)) for i in range(sites-1) for s in (0,1)),*(mono(((1,2*(i+1)+s),(0,2*i+s)),F(-1,3)) for i in range(sites-1) for s in (0,1)))
  m={'denominator':2,'terms':[{'powers':[2 if i==j else 0 for i in range(sites)],'coefficient':1} for j in range(sites)]} if metric=='D' else {'denominator':1,'terms':[{'powers':[0]*sites,'coefficient':1}]}
  c={'modes':modes,'particles':sites,'hamiltonian':encode(h),'polynomial_metric':m};start=time.monotonic()
  with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions allowed')):w,k,cost=JointPolynomial(c).compile(F(2))
  rows.append({'sites':sites,'metric':metric,'onsite_U':'4','hopping_t':'1/3','gamma':'2','seconds':time.monotonic()-start,'compilation':cost})
  if metric=='D':assert cost['numerator_degree']<=4 and cost['projector_free_zero_metric']
  else:assert cost['numerator_degree']==sites
out=Path(__file__).resolve().parents[1]/'projector_scaling_diagnostic.json';out.write_text(json.dumps({'scope':'Symbolic compilation only; no positivity or energy certificate at larger sizes.','rows':rows},indent=2)+'\n');print(out)
