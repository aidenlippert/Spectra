"""Export and replay V8 exact certificates without changing historical results."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json
from .v7_headroom import model,TOL
from .v7_certificate import Generator,Piece,check_certificate
from .v7_adaptive_taylor import adaptive_taylor
from .v7_verify import packed
from .v8_integer_taylor import fraction_free_taylor
from .v8_fractional_cover import adaptive_cover,check_evolution_cover,check_fractional_cover
ROOT=Path(__file__).resolve().parents[1]

def run():
 start=perf_counter();archive=[];failures=[];exact_pairs=0
 frozen=json.loads((ROOT/'results/v8/cost_replay.json').read_text())
 expected={(x['arm'],x['n'],x['family'],x['gamma'],x['horizon']):x for x in frozen['rows'] if x['repeat']==0}
 for n in (3,4,6):
  for family in ('xxz','mixed'):
   for gamma in (F(0),F(1,5),F(2)):
    for T in (F(1,5),F(1,2)):
     h,o=model(n,family);fraction_piece=None
     for arm in ('fraction','integer','partition','overlap'):
      want=expected[(arm,n,family,str(gamma),str(T))]
      try:
       g=Generator(h,gamma,n,512)
       if arm=='fraction':p,w,_=adaptive_taylor(g,o,T,TOL)
       elif arm=='integer':p,w,_=fraction_free_taylor(g,o,T,TOL)
       else:p,w,_=adaptive_cover(g,o,T,TOL,enabled=arm=='overlap')
      except ValueError as exc:
       if want['status']!='refused' or want['reason']!=str(exc):raise AssertionError('refusal drift')
       failures.append(dict(arm=arm,n=n,family=family,gamma=str(gamma),time=str(T),reason=str(exc)))
       continue
      if want['status']!='certified' or want['bound']!=w['claimed_bound']:raise AssertionError('result drift')
      if arm=='fraction':fraction_piece=p
      if arm=='integer':
       if p!=fraction_piece:raise AssertionError('integer coefficient mismatch')
       exact_pairs+=1
      archive.append(packed(n,h,gamma,o,T,p,w,dict(arm=arm,family=family)))
 archive_path=ROOT/'results/v8/certificate_archive.json'
 archive_path.write_text(json.dumps(archive,separators=(',',':'))+'\n')
 replays=[]
 for record in json.loads(archive_path.read_text()):
  g=Generator({p:F(c) for p,c in record['h'].items()},F(record['gamma']),record['n'],512)
  o={p:F(c) for p,c in record['initial'].items()}
  pieces=[Piece(F(p['duration']),tuple({label:F(c) for label,c in op.items()} for op in p['coefficients'])) for p in record['pieces']]
  check=check_certificate if record['tag']['arm'] in ('fraction','integer') else check_evolution_cover
  answer=check(g,o,pieces,record['witness'],F(record['tolerance']),expected_time=F(record['time']))
  if answer['status']!='certified':raise AssertionError(answer)
  replays.append(dict(tag=record['tag'],n=record['n'],gamma=record['gamma'],time=record['time'],bound=answer['bound']))
 norm_count=0
 probe=json.loads((ROOT/'results/v8/fractional_probe.json').read_text())
 for row in probe['rows']:
  if row['status']=='certified_norm':
   ans=check_fractional_cover(row['residual'],row['cover'])
   if ans['status']!='certified' or F(ans['bound'])*F(row['weight'])!=F(row['integrated_bound']):raise AssertionError('norm replay failed')
   norm_count+=1
 for source,expected_hash in {**frozen['source_sha256'],**probe['source_sha256']}.items():
  if sha256((ROOT/source).read_bytes()).hexdigest()!=expected_hash:raise AssertionError('source drift: '+source)
 v6=json.loads((ROOT/'results/v6/verification_receipt.json').read_text());v5=json.loads((ROOT/'results/v5/verification_receipt.json').read_text());v7=json.loads((ROOT/'results/v7/verification_receipt.json').read_text())
 historical={**v6['source_hashes'],**v6['data_hashes'],'results/v6/headroom_results.json':v6['headroom_result_sha256'],**v5['source_hashes'],**v5['result_hashes'],**v7['source_hashes']}
 checks={p:sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in historical.items()}
 if not all(checks.values()):raise AssertionError('historical drift')
 if sha256((ROOT/'results/v7/certificate_archive.json').read_bytes()).hexdigest()!=v7['archive_sha256']:raise AssertionError('V7 archive drift')
 files=list((ROOT/'experiments').glob('v8*.py'))+list((ROOT/'tests').glob('test_v8*.py'))
 report=dict(status='verified',evolution_certificates_replayed=len(replays),norm_certificates_replayed=norm_count,
             exact_fraction_integer_piece_pairs=exact_pairs,refusals_replayed=failures,receipts=replays,
             historical_hash_checks=checks,v7_archive_unchanged=True,archive_sha256=sha256(archive_path.read_bytes()).hexdigest(),
             source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in files},
             elapsed_seconds=perf_counter()-start,compounding_capability='not_established',acquisition='not_run',physical_validation='not_performed')
 (ROOT/'results/v8/verification_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
 print({k:report[k] for k in ('status','evolution_certificates_replayed','norm_certificates_replayed','exact_fraction_integer_piece_pairs','elapsed_seconds')})
if __name__=='__main__':run()
