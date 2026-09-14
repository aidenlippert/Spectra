"""Whole approximation construction diagnostic, no learning or heldout use."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import hashlib,json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,check_certificate
from .v8_integer_taylor import fraction_free_taylor
from .v9_galerkin import adaptive_galerkin
from .v9_polynomial import krylov_basis,collocation_piece,projected_taylor_piece,best_witness
from .v7_verify import packed
ROOT=Path(__file__).resolve().parents[1]

def attempt(n,family,gamma,T,arm):
 start=perf_counter();trials=[];g=None
 try:
  h,o=model(n,family);g=Generator(h,gamma,n,512)
  if arm=='integer':
   p,w,c=fraction_free_taylor(g,o,T,TOL)
   receipt=check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=T)
   if receipt['status']!='certified':raise ValueError('checker: '+str(receipt))
  elif arm in ('galerkin','bfs'):
   c=adaptive_galerkin(g,o,T,TOL,growth='bfs' if arm=='bfs' else 'residual')
   if c['status']!='certified':return dict(status='refused',reason=c['reason'],cost=c,complete_seconds=perf_counter()-start),None
   p,w,receipt=c.pop('piece'),c.pop('certificate'),c['checker']
  else:
   accepted=None
   for dimension in (8,16):
    bstart=perf_counter();basis=krylov_basis(g,o,dimension);basis_seconds=perf_counter()-bstart
    degrees=(4,6,8,10,12,16,20,24) if arm=='krylov_taylor' else (4,6,8,10,12)
    for degree in degrees:
     t0=perf_counter();construct=projected_taylor_piece if arm=='krylov_taylor' else collocation_piece
     p,pc=construct(basis,o,T,degree)
     w,wc=best_witness(Generator(h,gamma,n,512),o,p,TOL)
     row=dict(dimension=dimension,degree=degree,proposed_bound=w['claimed_bound'],basis_seconds=basis_seconds if degree==degrees[0] else 0.,basis_cost=basis['cost'] if degree==degrees[0] else None,polynomial_work=pc,witness_work=wc,attempt_seconds=perf_counter()-t0)
     trials.append(row)
     if F(w['claimed_bound'])<=TOL:
      ct=perf_counter();receipt=check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=T)
      row['checking_seconds']=perf_counter()-ct;row['checker']=receipt
      if receipt['status']!='certified':raise ValueError('checker: '+str(receipt))
      accepted=(p,w);break
    if accepted:break
   if not accepted:return dict(status='refused',reason='bounded polynomial schedule exhausted',trials=trials,generator=dict(g.cost),complete_seconds=perf_counter()-start),None
   c=dict(trials=trials,generator=dict(g.cost))
  artifact=packed(n,h,gamma,o,T,p,w,dict(arm=arm,family=family))
  return dict(status='certified',bound=w['claimed_bound'],order=len(p.coefficients)-1,entries=sum(map(len,p.coefficients)),cost=c,checker=receipt,complete_seconds=perf_counter()-start),artifact
 except (ValueError,RuntimeError) as exc:
  return dict(status='refused',reason=str(exc),trials=trials,generator=dict(g.cost) if g else {},complete_seconds=perf_counter()-start),None

def run():
 rows=[];archive=[]
 for n in (3,4,6):
  for family in ('xxz','mixed'):
   for gamma in (F(0),F(1,5),F(2)):
    for T in (F(1,5),F(1,2)):
     for arm in ('integer','galerkin','bfs','krylov_taylor','collocation'):
      row,artifact=attempt(n,family,gamma,T,arm);row.update(n=n,family=family,gamma=str(gamma),time=str(T),arm=arm);rows.append(row)
      if artifact:archive.append(artifact)
      print(n,family,gamma,T,arm,row['status'],round(row['complete_seconds'],3),flush=True)
 files=['research/v9/HEADROOM_PROTOCOL.md','experiments/v9_headroom.py','experiments/v9_polynomial.py','experiments/v9_galerkin.py','experiments/v8_integer_taylor.py','experiments/v7_certificate.py']
 out=dict(schema='v9-development-1',rows=rows,source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},scope='development only; supplied conventional methods, no acquired method')
 (ROOT/'results/v9/headroom.json').write_text(json.dumps(out,indent=2)+'\n')
 (ROOT/'results/v9/certificate_archive.json').write_text(json.dumps(archive,separators=(',',':'))+'\n')
 print('complete',len(rows),'accepted',len(archive),flush=True)
if __name__=='__main__':run()
