from fractions import Fraction as F
from time import perf_counter
from random import Random
from hashlib import sha256
from statistics import median
import json
from .v7_certificate import Generator
from .v7_quadratic_certificate import adaptive_with_square,check
from .v7_headroom import ROOT,model,TOL

def run():
 rng=Random(271828);rows=[]
 for n in (3,4,6):
  for family in ('xxz','mixed'):
   for gamma in (F(0),F(1,5),F(2)):
    for T in (F(1,5),F(1,2)):
     case=[]
     for repeat in range(7):
      arms=[1,40];rng.shuffle(arms)
      for cap in arms:
       start=perf_counter();h,o=model(n,family);g=Generator(h,gamma,n,max_terms=512)
       r=dict(n=n,family=family,gamma=str(gamma),T=str(T),repeat=repeat,arm='baseline' if cap==1 else 'quadratic')
       try:
        p,w,c=adaptive_with_square(g,o,T,TOL,term_cap=cap)
        chk=check(Generator(h,gamma,n,max_terms=512),o,[p],w,TOL,expected_time=T)
        r.update(status=chk['status'],order=len(p.coefficients)-1,bound=chk.get('bound'),quadratic_tests=c['quadratic_tests'],checking_cost=chk['cost'])
        if chk['status']!='certified':raise AssertionError(chk)
       except ValueError as exc:r.update(status='no_certificate_within_budget',reason=str(exc))
       r['elapsed']=perf_counter()-start;rows.append(r);case.append(r)
     successful=[r for r in case if r['status']=='certified']
     if successful:
      b=median([r['elapsed'] for r in case if r['arm']=='baseline']);q=median([r['elapsed'] for r in case if r['arm']=='quadratic'])
      print(n,family,str(gamma),str(T),'ratio',round(b/q,3),flush=True)
 result=dict(rows=rows,protocol_sha256=sha256((ROOT/'research/v7/QUADRATIC_PROTOCOL.md').read_bytes()).hexdigest(),checker_sha256=sha256((ROOT/'experiments/v7_quadratic_certificate.py').read_bytes()).hexdigest(),autonomous_acquisition='not_run')
 (ROOT/'results/v7/quadratic_cost.json').write_text(json.dumps(result,indent=2));return result
if __name__=='__main__':run()
