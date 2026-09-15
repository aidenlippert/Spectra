"""Development-only norm probe; no full-cost headroom or method-learning claim."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import hashlib,json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,norm_witness
from .v7_adaptive_taylor import adaptive_taylor
from .v8_fractional_cover import propose_fractional_cover,check_fractional_cover

ROOT=Path(__file__).resolve().parents[1]

def run():
 rows=[]
 for n in (3,4):
  for family in ('xxz','mixed'):
   for gamma in (F(0),F(1,5),F(2)):
    for T in (F(1,5),F(1,2)):
     h,initial=model(n,family)
     p,w,_=adaptive_taylor(Generator(h,gamma,n,max_terms=512),initial,T,TOL)
     for degree in (len(p.coefficients)-2,len(p.coefficients)-1):
      derivative=Generator(h,gamma,n,max_terms=512).apply(p.coefficients[degree])
      factor=T**(degree+1)/F(degree+1)
      start=perf_counter();partitions={}
      for mode in ('l1','firstfit','weighted'):
       groups,cost=norm_witness(derivative,mode)
       partitions[mode]=dict(bound=str(sum((F(g['upper']) for g in groups),F(0))*factor),cost=cost)
      baseline_seconds=perf_counter()-start
      start=perf_counter()
      try:
       cover,cost=propose_fractional_cover(derivative,iterations=16)
       receipt=check_fractional_cover(derivative,cover)
       if receipt['status']!='certified':raise ValueError(receipt['reason'])
       bound=F(receipt['bound'])*factor
       row=dict(status='certified_norm',cover=cover,cost=cost,checker=receipt,integrated_bound=str(bound),
                newly_meets_tolerance=bound<=TOL and min(F(v['bound']) for v in partitions.values())>TOL)
      except ValueError as e:row=dict(status='refused',reason=str(e))
      row.update(n=n,family=family,gamma=str(gamma),horizon=str(T),order=degree,terms=len(derivative),
                 comparison_partitions=partitions,partition_seconds=baseline_seconds,
                 cover_total_seconds=perf_counter()-start,residual={p:str(c) for p,c in derivative.items()},weight=str(factor))
      rows.append(row)
     print(n,family,gamma,T,'done',flush=True)
 out=dict(schema='v8-fractional-probe-2',scope='existing development residuals only; no complete calculation or acquisition claim',
          tolerance=str(TOL),rows=rows,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'experiments/v8_fractional_cover.py')})
 (ROOT/'results/v8/fractional_probe.json').write_text(json.dumps(out,indent=2)+'\n')
 print('rows',len(rows),'new_tolerance',sum(r.get('newly_meets_tolerance',False) for r in rows),flush=True)
if __name__=='__main__':run()
