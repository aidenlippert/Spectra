"""Costed conventional development calculations. No learner or held-out access."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from hashlib import sha256
import json
from .v7_certificate import Generator, Piece, derive_certificate, check_certificate
from .v7_reducers import full_taylor, projected_taylor, bfs_basis, residual_basis, arnoldi
ROOT=Path(__file__).resolve().parents[1]
TOL=F(1,1000)

def model(n,family):
 h={}
 def put(sites,value):
  p=['I']*n
  for i,letter in sites:p[i]=letter
  word=''.join(p);h[word]=h.get(word,F(0))+value
 for i in range(n):
  put([(i,'Z')],F(1,2))
  if family=='mixed':put([(i,'X')],F((-1)**i,5))
 for i in range(n-1):
  put([(i,'X'),(i+1,'X')],F(3+i%2,5))
  put([(i,'Z'),(i+1,'Z')],F(1,5))
  if family=='mixed':put([(i,'Y'),(i+1,'Y')],F(1,3))
 initial={''.join('Z' if i==n//2 else 'I' for i in range(n)):F(1)}
 return h,initial

def settings(method):
 if method=='taylor':return [(p,None) for p in (4,8,12,16,24)]
 if method=='bfs':return [(p,d) for p,d in ((8,2),(12,4),(16,6),(24,10))]
 if method=='residual':return [(p,b) for p,b in ((8,2),(12,8),(16,32),(24,128))]
 if method=='arnoldi':return [(p,k) for p,k in ((8,4),(12,8),(16,12),(24,20))]
 raise ValueError('method')

def calculation(n,family,gamma,T,method,deadline=20):
 began=perf_counter();h,initial=model(n,family);attempts=[];phases=dict(construction_solve=0.,witness=0.,checking=0.)
 certified=False;best=None
 for order,extra in settings(method):
  gen=Generator(h,gamma,n,max_terms=512);start=perf_counter()
  try:
   if method=='taylor':piece,work=full_taylor(gen,initial,T,order)
   elif method=='bfs':
    basis,_=bfs_basis(gen,initial,extra,cap=512)
    piece,work=projected_taylor(gen,initial,T,order,basis)
   elif method=='residual':piece,work=residual_basis(gen,initial,T,order,cap=512,batch=extra)
   else:piece,work=arnoldi(gen,initial,T,order,extra)
   phases['construction_solve']+=perf_counter()-start
   for integration in ('power','bernstein'):
    for grouping in ('l1','firstfit','weighted'):
     if perf_counter()-began>deadline:raise RuntimeError('calculation time budget')
     start=perf_counter();wgen=Generator(h,gamma,n,max_terms=512)
     witness=derive_certificate(wgen,initial,[piece],grouping,integration)
     phases['witness']+=perf_counter()-start
     start=perf_counter();cgen=Generator(h,gamma,n,max_terms=512)
     result=check_certificate(cgen,initial,[piece],witness,TOL,expected_time=T)
     phases['checking']+=perf_counter()-start
     row=dict(order=order,extra=extra,grouping=grouping,integration_basis=integration,
              status=result['status'],bound=result.get('bound'),construction_work=work,
              witness_work=witness['construction_cost'],witness_generator=witness['generator_cost'],checking=result,
              polynomial_entries=sum(map(len,piece.coefficients)),peak_support=max(map(len,piece.coefficients)))
     attempts.append(row)
     if result['status']=='rejected':raise AssertionError(result['reason'])
     if result['status']=='certified':certified=True;best=row;break
    if certified:break
   if certified:break
  except (ValueError,RuntimeError) as exc:
   # Failed attempt work and elapsed time remain in the accumulated account.
   attempts.append(dict(order=order,extra=extra,status='refused',reason=str(exc),generator_work=dict(gen.cost)))
   if 'time budget' in str(exc):break
  if perf_counter()-began>deadline:break
 return dict(n=n,family=family,gamma=str(gamma),T=str(T),method=method,tolerance=str(TOL),
             status='certified' if certified else 'no_certificate_within_budget',
             elapsed=perf_counter()-began,phases=phases,attempts=attempts,accepted=best)

def run(smoke=False):
 protocol=ROOT/'research/v7/HEADROOM_PROTOCOL.md';rows=[]
 sizes=(3,) if smoke else (3,4,6)
 families=('xxz',) if smoke else ('xxz','mixed')
 for n in sizes:
  for family in families:
   for gamma in (F(0),F(1,5),F(2)):
    for T in ((F(1,5),) if smoke else (F(1,5),F(1,2))):
     for method in ('taylor','bfs','residual','arnoldi'):
      row=calculation(n,family,gamma,T,method);rows.append(row)
      print(n,family,str(gamma),str(T),method,row['status'],round(row['elapsed'],4),flush=True)
 result=dict(protocol_sha256=sha256(protocol.read_bytes()).hexdigest(),checker_sha256=sha256((ROOT/'experiments/v7_certificate.py').read_bytes()).hexdigest(),
             smoke=smoke,rows=rows,headroom_gate='not_established_conventional_baselines_only',autonomous_acquisition='not_run')
 out=ROOT/('results/v7/headroom_smoke.json' if smoke else 'results/v7/headroom_development.json')
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2))
 return result
if __name__=='__main__':
 import sys
 run('--smoke' in sys.argv)
