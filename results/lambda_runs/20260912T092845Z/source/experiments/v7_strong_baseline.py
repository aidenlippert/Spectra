from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from hashlib import sha256
import json
from .v7_certificate import Generator, check_certificate
from .v7_adaptive_taylor import adaptive_taylor
from .v7_headroom import model, TOL, ROOT

def run():
 rows=[]
 for n in (3,4,6):
  for family in ('xxz','mixed'):
   for gamma in (F(0),F(1,5),F(2)):
    for T in (F(1,5),F(1,2)):
     start=perf_counter();h,o=model(n,family);g=Generator(h,gamma,n,max_terms=512)
     row=dict(n=n,family=family,gamma=str(gamma),T=str(T),method='adaptive_taylor',tolerance=str(TOL))
     try:
      p,w,c=adaptive_taylor(g,o,T,TOL);solve_end=perf_counter()
      checker=Generator(h,gamma,n,max_terms=512)
      receipt=check_certificate(checker,o,[p],w,TOL,expected_time=T)
      row.update(status=receipt['status'],bound=receipt.get('bound'),construction_solve_proposal_elapsed=solve_end-start,
                 checking_elapsed=perf_counter()-solve_end,construction=c,checking=receipt,order=len(p.coefficients)-1,
                 polynomial_entries=sum(map(len,p.coefficients)),witness=w)
     except ValueError as exc:
      row.update(status='no_certificate_within_budget',reason=str(exc),construction_generator=dict(g.cost))
     row['elapsed']=perf_counter()-start;rows.append(row)
     print(n,family,str(gamma),str(T),row['status'],round(row['elapsed'],4),flush=True)
 result=dict(rows=rows,headroom_gate='not_established_conventional_improvements_incorporated',
             source_hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ('experiments/v7_certificate.py','experiments/v7_adaptive_taylor.py','research/v7/BASELINE_AMENDMENT.md')})
 (ROOT/'results/v7/strong_baseline.json').write_text(json.dumps(result,indent=2))
 return result
if __name__=='__main__':run()
