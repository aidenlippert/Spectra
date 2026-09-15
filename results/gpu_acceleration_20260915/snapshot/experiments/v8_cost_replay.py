"""Paired development-only complete-cost comparisons, never an acquisition run."""
from fractions import Fraction as F
from pathlib import Path
from random import Random
from time import perf_counter
import json,hashlib
from .v7_headroom import model,TOL
from .v7_certificate import Generator,check_certificate
from .v7_adaptive_taylor import adaptive_taylor
from .v8_integer_taylor import fraction_free_taylor
from .v8_fractional_cover import adaptive_cover,check_evolution_cover
ROOT=Path(__file__).resolve().parents[1]

def attempt(n,family,gamma,T,arm):
 start=perf_counter();gen=None
 try:
  h,o=model(n,family);gen=Generator(h,gamma,n,512)
  if arm=='fraction':p,w,c=adaptive_taylor(gen,o,T,TOL)
  elif arm=='integer':p,w,c=fraction_free_taylor(gen,o,T,TOL)
  else:p,w,c=adaptive_cover(gen,o,T,TOL,enabled=(arm=='overlap'))
  construct=perf_counter()-start;t0=perf_counter();cg=Generator(h,gamma,n,512)
  checker=check_certificate if arm in ('fraction','integer') else check_evolution_cover
  receipt=checker(cg,o,[p],w,TOL,expected_time=T)
  checked=perf_counter()-t0
  return dict(status=receipt['status'],bound=receipt.get('bound'),reason=receipt.get('reason'),order=len(p.coefficients)-1,
              construction_seconds=construct,checking_seconds=checked,complete_seconds=perf_counter()-start,
              cost=c,checker=receipt,entries=sum(map(len,p.coefficients)))
 except ValueError as exc:
  return dict(status='refused',reason=str(exc),complete_seconds=perf_counter()-start,generator=dict(gen.cost) if gen else {})

def run():
 rng=Random(8031);rows=[]
 for repeat in range(5):
  for n in (3,4,6):
   for family in ('xxz','mixed'):
    for gamma in (F(0),F(1,5),F(2)):
     for T in (F(1,5),F(1,2)):
      arms=['fraction','integer','partition','overlap'];rng.shuffle(arms)
      for arm in arms:
       row=attempt(n,family,gamma,T,arm);row.update(repeat=repeat,n=n,family=family,gamma=str(gamma),horizon=str(T),arm=arm);rows.append(row)
  print('repeat',repeat,'complete',flush=True)
 files=['research/v8/COST_PROTOCOL.md','experiments/v8_cost_replay.py','experiments/v8_integer_taylor.py','experiments/v8_fractional_cover.py','experiments/v7_certificate.py','experiments/v7_adaptive_taylor.py','experiments/v7_headroom.py','experiments/certificates.py','experiments/pauli.py']
 out=dict(schema='v8-cost-replay-1',rows=rows,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},scope='development, no learned method and no heldout evaluation')
 (ROOT/'results/v8/cost_replay.json').write_text(json.dumps(out,indent=2)+'\n')
 print('rows',len(rows),flush=True)
if __name__=='__main__':run()
