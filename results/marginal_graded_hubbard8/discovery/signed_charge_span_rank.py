"""Exact rank and dependencies before enlarging the signed-charge telescope span."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from math import prod
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from experiments.marginal_quadratic_charge_telescope import LABELS as Q,five_site_value as qvalue
from experiments.marginal_charge_square_pairs import LABELS as S,five_site_value as svalue
from experiments.marginal_charge_indicator_telescope import LABELS as I,five_site_value as ivalue
source=ROOT/'results/marginal_graded_hubbard8/full_charge_indicators/W_zero/polished/range_two_family_limit_certificate.json'
c=json.loads(source.read_text());shapes=[{int(s):F(v) for s,v in d.items()} for d in c['diagonal_shapes']]
exponents=[e for e in product(range(3),repeat=5) if sum(e)%2==0 and e<e[::-1] and sum(e)!=2 and 1 in e]
assert len(exponents)==36
pivots={}
for state in range(1024):
 q=[((state>>(2*i))&3).bit_count()-1 for i in range(5)]
 row=[d.get(state,F(0)) for d in shapes]+[fn(state,{key:F(1)}) for names,fn in [(Q,qvalue),(S,svalue),(I,ivalue)] for key in names]+[F(prod(x**power for x,power in zip(q,e))-prod(x**power for x,power in zip(q,e[::-1]))) for e in exponents]
 for j in sorted(pivots):
  a=row[j]
  if a:row=[v-a*b for v,b in zip(row,pivots[j])]
 j=next((j for j,v in enumerate(row) if v),None)
 if j is not None:
  a=row[j];pivots[j]=[v/a for v in row]
files={Path(__file__).resolve(),source}
for module in tuple(sys.modules.values()):
 path=getattr(module,'__file__',None)
 if path and str(Path(path).resolve()).startswith(str(ROOT/'experiments')+'/'):files.add(Path(path).resolve())
result={'accepted':True,'physical_five_site_determinants':1024,'baseline_columns':25,'mixed_columns':36,'total_columns':61,'rank':len(pivots),'independent_mixed_indices':[i for i in range(36) if i+25 in pivots],'exponents':exponents,'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},'scope':'Exact rank of the full five-site diagonal operator columns; no energy or family-limit acceptance. Reflection-odd Y has no constant component, so Y_left-Y_right has the same rank.'}
(ROOT/'results/marginal_graded_hubbard8/full_signed_charge/span_rank.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','exponents')}))
