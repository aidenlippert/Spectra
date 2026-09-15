import json,itertools,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
ROOT=Path('/Users/aidenlippert/Documents/Spectra'); OUT=ROOT/'results/wave2_20260913/geometry'; OUT.mkdir(parents=True,exist_ok=True)
def jw(word,m,inds):
 n=1<<len(inds); M=np.zeros((n,n)); pos={x:i for i,x in enumerate(inds)}
 for bits in range(n):
  col=bits; amp=1
  for typ,p in word[::-1]:
   if p not in pos: return M
   j=pos[p]; occ=(col>>j)&1; sign=(-1)**(sum((col>>k)&1 for k in range(j)))
   if typ==0:
    if not occ: amp=0; break
    amp*=sign; col^=1<<j
   else:
    if occ: amp=0; break
    amp*=sign; col^=1<<j
  if amp: M[col,bits]+=amp
 return M
def run(path,clusters,label):
 d=json.loads(path.read_text()); m=d['modes']; blocks=[[] for _ in clusters]; residual=Fraction(0); assigned=0
 for row in d['hamiltonian']:
  w=tuple(tuple(x) for x in row['word']); c=Fraction(row['coefficient']); supp=set(x[1] for x in w)
  hit=next((i for i,C in enumerate(clusters) if supp<=set(C)),None)
  if hit is None: residual+=abs(c)
  else: blocks[hit].append((w,c)); assigned+=1
 lows=[]; costs=[]
 for C,br in zip(clusters,blocks):
  K=len(C); A=np.zeros((1<<K,1<<K),dtype=object)
  for w,c in br: A += object_matrix(jw(w,m,C),c)
  herm=all(A[i,j]==A[j,i] for i in range(1<<K) for j in range(1<<K))
  if not herm: raise ValueError('non-Hermitian local block')
  # Gershgorin lower bound is exact rational: min diagonal - row absolute offdiag.
  gl=[]
  for i in range(1<<K):
   rad=sum(abs(A[i,j]) for j in range(1<<K) if j!=i); gl.append(A[i,i]-rad)
  lo=min(gl); lows.append(lo); costs.append({'cluster':C,'terms':len(br),'dimension':1<<K,'gershgorin_lower_exact':str(lo)})
 lb=sum(lows)-residual
 return {'variant':label,'fixture':str(path.relative_to(ROOT)),'fixture_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'clusters':clusters,'max_cluster_width':max(map(len,clusters)),'assigned_terms':assigned,'total_terms':len(d['hamiltonian']),'omitted_residual_l1_exact':str(residual),'block_cost':costs,'all_blocks_hermitian':True,'sum_block_lower_exact':str(sum(lows)),'global_lower_exact':str(lb),'global_lower':float(lb),'assignment':'first fitting cluster, weight 1; fitting clusters receive zero weight; nonfitting words charged L1'}
def object_matrix(M,c): return np.vectorize(lambda x: c*int(round(x)),otypes=[object])(M)
res=[]
for i in (4,6):
 p=ROOT/f'results/certificate_scaling/active_space_ladder/h{i}/fixture.json'; m=2*i
 res.append(run(p,[list(range(0,4)),list(range(4,8)),list(range(8,m))],f'H{i}_disjoint4'))
 windows=[list(range(0,6)),list(range(4,min(m,10))),list(range(6,m))] if m>=12 else [list(range(0,6)),list(range(2,8))]
 res.append(run(p,windows,f'H{i}_overlap6'))
(OUT/'overlap_cluster_lower_receipt.json').write_text(json.dumps(res,indent=2)+'\n'); print(json.dumps(res,indent=2))
