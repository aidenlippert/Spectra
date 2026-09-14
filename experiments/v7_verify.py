"""Regenerate exact accepted artifacts and independently replay their checks."""
from fractions import Fraction as F
from hashlib import sha256
from time import perf_counter
import json
from .v7_certificate import Generator,Piece,check_certificate,derive_certificate
from .v7_adaptive_taylor import adaptive_taylor
from .v7_reducers import full_taylor,projected_taylor,bfs_basis,residual_basis,arnoldi
from .v7_headroom import ROOT,model,TOL
from .v7_measured_gate import run as measured

def packed(n,h,gamma,initial,T,piece,witness,tag):
 return dict(tag=tag,n=n,h={p:str(v) for p,v in h.items()},gamma=str(gamma),initial={p:str(v) for p,v in initial.items()},time=str(T),tolerance=str(TOL),
             pieces=[dict(duration=str(piece.duration),coefficients=[{p:str(v) for p,v in c.items()} for c in piece.coefficients])],witness=witness)

def replay(record):
 g=Generator({p:F(v) for p,v in record['h'].items()},F(record['gamma']),record['n'],max_terms=512)
 initial={p:F(v) for p,v in record['initial'].items()}
 pieces=[Piece(F(r['duration']),tuple({p:F(v) for p,v in c.items()} for c in r['coefficients'])) for r in record['pieces']]
 return check_certificate(g,initial,pieces,record['witness'],F(record['tolerance']),expected_time=F(record['time']))

def run():
 start=perf_counter();archive=[];failure_replays=0
 strong=json.loads((ROOT/'results/v7/strong_baseline.json').read_text())
 for r in strong['rows']:
  n=r['n'];h,o=model(n,r['family']);gamma=F(r['gamma']);T=F(r['T'])
  try:p,w,_=adaptive_taylor(Generator(h,gamma,n,max_terms=512),o,T,TOL)
  except ValueError:
   if r['status']=='certified':raise AssertionError('saved success no longer reproducible')
   failure_replays+=1;continue
  if r['status']!='certified' or w['claimed_bound']!=r['bound']:raise AssertionError('strong baseline drift')
  archive.append(packed(n,h,gamma,o,T,p,w,dict(arm='adaptive_taylor',family=r['family'])))
 dev=json.loads((ROOT/'results/v7/headroom_development.json').read_text())
 for r in dev['rows']:
  if r['status']!='certified':continue
  n=r['n'];h,o=model(n,r['family']);gamma=F(r['gamma']);T=F(r['T']);g=Generator(h,gamma,n,max_terms=512)
  a=r['accepted'];order=a['order'];extra=a['extra']
  if r['method']=='taylor':p,_=full_taylor(g,o,T,order)
  elif r['method']=='bfs':
   b,_=bfs_basis(g,o,extra,cap=512);p,_=projected_taylor(g,o,T,order,b)
  elif r['method']=='residual':p,_=residual_basis(g,o,T,order,cap=512,batch=extra)
  else:p,_=arnoldi(g,o,T,order,extra)
  w=derive_certificate(Generator(h,gamma,n,max_terms=512),o,[p],a['grouping'],a['integration_basis'])
  if w['claimed_bound']!=a['bound']:raise AssertionError('development certificate drift')
  archive.append(packed(n,h,gamma,o,T,p,w,dict(arm=r['method'],family=r['family'])))
 archive_path=ROOT/'results/v7/certificate_archive.json';archive_path.write_text(json.dumps(archive,separators=(',',':')))
 # Reparse exported bytes before invoking the independent checker.
 receipts=[]
 for r in json.loads(archive_path.read_text()):
  ans=replay(r)
  if ans['status']!='certified':raise AssertionError(ans)
  receipts.append(dict(tag=r['tag'],n=r['n'],gamma=r['gamma'],time=r['time'],bound=ans['bound']))
 mr=measured();saved=json.loads((ROOT/'results/v7/prediction_requests.json').read_text())
 if mr!=saved:raise AssertionError('measured gate drift')
 # Previous versions are references, never regenerated or rewritten here.
 previous=json.loads((ROOT/'results/v6/verification_receipt.json').read_text())
 hashes={**previous['source_hashes'],**previous['data_hashes'],
         'results/v6/headroom_results.json':previous['headroom_result_sha256']}
 v5=json.loads((ROOT/'results/v5/verification_receipt.json').read_text())
 hashes.update(v5['source_hashes']);hashes.update(v5['result_hashes'])
 preserved={p:sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items()}
 if not all(preserved.values()):raise AssertionError('historical hash changed')
 sources=list((ROOT/'experiments').glob('v7*.py'))+[ROOT/'experiments/pauli.py',ROOT/'experiments/certificates.py']
 report=dict(status='verified',accepted_certificates_replayed=len(receipts),strong_baseline_failure_replays=failure_replays,
             measured_witness_replayed=True,archive_sha256=sha256(archive_path.read_bytes()).hexdigest(),archive_bytes=archive_path.stat().st_size,
             source_hashes={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
             elapsed_seconds=perf_counter()-start,receipts=receipts,historical_hash_checks=preserved,
             compounding_capability='not_established',autonomous_acquisition='not_run',v4_calibration_dependency='excluded from V5',physical_validation='not_performed')
 (ROOT/'results/v7/verification_receipt.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:report[k] for k in ('status','accepted_certificates_replayed','strong_baseline_failure_replays','archive_bytes','elapsed_seconds')},indent=2));return report
if __name__=='__main__':run()
