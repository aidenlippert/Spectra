"""Exact geometry-aware lower bound: diagonal local blocks plus residual L1.

The retained blocks are diagonal in occupation and partitioned by spatial
orbital (spin pairs). Their minima are enumerated in the requested particle
sector. Every omitted normal-ordered CAR coefficient is charged by |c|,
which is a rigorous operator-norm bound (each monomial has norm <= 1).
"""
import json,itertools
from pathlib import Path
from fractions import Fraction
ROOT=Path('/Users/aidenlippert/Documents/Spectra')
FIX=[ROOT/f'results/certificate_scaling/active_space_ladder/h{i}/fixture.json' for i in (4,6)]
OUT=ROOT/'results/wave2_20260913/geometry'; OUT.mkdir(parents=True,exist_ok=True)
def diagonal(w):
    if len(w)==2: return w[0][0]==1 and w[1][0]==0 and w[0][1]==w[1][1]
    if len(w)==4 and [x[0] for x in w]==[1,1,0,0]:
        return sorted((w[0][1],w[1][1]))==sorted((w[2][1],w[3][1]))
    return False
def occ(w,b):
    # evaluate normal ordered diagonal monomial on occupation bitstring
    if len(w)==2:return b[w[0][1]]
    p,q=w[0][1],w[1][1]
    return b[p]*b[q] if p!=q else 0
def run(path):
 d=json.loads(path.read_text()); m=d['modes']; N=d['particles']; rows=d['hamiltonian']
 # spatial clusters are spin pairs; for H4/H6 these are contiguous alpha,beta pairs
 clusters=[list(range(2*i,min(2*i+2,m))) for i in range(m//2)]
 block_rows=[[] for _ in clusters]; residual=Fraction(0); retained=0
 for row in rows:
  w=tuple(tuple(x) for x in row['word']); c=Fraction(row['coefficient'])
  # local if diagonal and all occupied orbital indices in one spatial block
  inds=[x[1] for x in w]; bi={j:j//2 for j in inds}
  if diagonal(w) and len(set(bi.values()))==1:
   block_rows[next(iter(bi.values()))].append((w,c)); retained+=1
  else: residual+=abs(c)
 mins=[]; block_cost=[]
 for k,br in enumerate(block_rows):
  vals=[]
  for bits in itertools.product((0,1),repeat=m):
   if sum(bits)!=N: continue
   vals.append(sum(c*occ(w,bits) for w,c in br))
  mins.append(min(vals) if vals else Fraction(0)); block_cost.append({'cluster':clusters[k],'terms':len(br),'states_tested':len(vals)})
 lb=sum(mins)-residual
 return {'fixture':str(path.relative_to(ROOT)),'modes':m,'particles':N,'clusters':clusters,'block_cost':block_cost,'retained_diagonal_terms':retained,'total_terms':len(rows),'local_minima_exact':[str(x) for x in mins],'local_minimum_sum_exact':str(sum(mins)),'residual_l1_exact':str(residual),'lower_bound_exact':str(lb),'lower_bound':float(lb),'validity':'H = retained diagonal blocks + exact residual; each block minimum is sector-enumerated and ||CAR word||<=1'}
res=[run(p) for p in FIX]; (OUT/'local_collective_lower_receipt.json').write_text(json.dumps(res,indent=2)+'\n'); print(json.dumps(res,indent=2))
