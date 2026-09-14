from fractions import Fraction as F
from time import perf_counter
from hashlib import sha256
import json
from .v7_headroom import ROOT,model,TOL
from .v7_certificate import Generator,check_certificate
from .v7_pruned_taylor import pruned_taylor

def run():
 rows=[]
 for n in (3,4,6):
  for family in ('xxz','mixed'):
   for gamma in (F(0),F(1,5),F(2)):
    for T in (F(1,5),F(1,2)):
     for mode in ('magnitude','work'):
      started=perf_counter();h,o=model(n,family);g=Generator(h,gamma,n,max_terms=512)
      row=dict(n=n,family=family,gamma=str(gamma),T=str(T),mode=mode)
      try:
       p,w,c=pruned_taylor(g,o,T,TOL,mode=mode);proposed=perf_counter()
       checked=check_certificate(Generator(h,gamma,n,max_terms=512),o,[p],w,TOL,expected_time=T)
       row.update(status=checked['status'],bound=checked.get('bound'),construction=c,checking=checked,
                  construction_solve_witness_time=proposed-started,checking_time=perf_counter()-proposed,
                  order=len(p.coefficients)-1,polynomial_entries=sum(map(len,p.coefficients)))
      except ValueError as exc:row.update(status='no_certificate_within_budget',reason=str(exc),generator=dict(g.cost))
      row['elapsed']=perf_counter()-started;rows.append(row)
      print(n,family,str(gamma),str(T),mode,row['status'],round(row['elapsed'],4),flush=True)
 result=dict(rows=rows,protocol_sha256=sha256((ROOT/'research/v7/PRUNING_PROTOCOL.md').read_bytes()).hexdigest(),checker_sha256=sha256((ROOT/'experiments/v7_certificate.py').read_bytes()).hexdigest(),autonomous_acquisition='not_run')
 (ROOT/'results/v7/pruning_probe.json').write_text(json.dumps(result,indent=2));return result
if __name__=='__main__':run()
