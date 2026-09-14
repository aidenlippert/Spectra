"""Independent saved-certificate replay and source/preservation checks for V9."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from collections import Counter
import json
from .v7_verify import replay
ROOT=Path(__file__).resolve().parents[1]

def run():
 start=perf_counter();saved=json.loads((ROOT/'results/v9/headroom.json').read_text());archive=json.loads((ROOT/'results/v9/certificate_archive.json').read_text())
 expected={(r['arm'],r['n'],r['family'],r['gamma'],r['time']):r for r in saved['rows'] if r['status']=='certified'}
 if len(expected)!=len(archive):raise AssertionError('archive count')
 receipts=[]
 for record in archive:
  key=(record['tag']['arm'],record['n'],record['tag']['family'],record['gamma'],record['time'])
  if key not in expected:raise AssertionError('unknown/duplicate certificate')
  result=replay(record)
  if result['status']!='certified' or result['bound']!=expected.pop(key)['bound']:raise AssertionError(result)
  receipts.append(dict(arm=key[0],n=key[1],family=key[2],gamma=key[3],time=key[4],bound=result['bound']))
 if expected:raise AssertionError('missing certificate')
 for p,h in saved['source_sha256'].items():
  if sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise AssertionError('source drift '+p)
 v8=json.loads((ROOT/'results/v8/verification_receipt.json').read_text())
 for p,h in v8['source_sha256'].items():
  if sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise AssertionError('v8 source drift '+p)
 for version in ('v7','v8'):
  old=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
  if sha256((ROOT/f'results/{version}/certificate_archive.json').read_bytes()).hexdigest()!=old['archive_sha256']:raise AssertionError('prior archive drift')
 v7=json.loads((ROOT/'results/v7/verification_receipt.json').read_text());v6=json.loads((ROOT/'results/v6/verification_receipt.json').read_text());v5=json.loads((ROOT/'results/v5/verification_receipt.json').read_text())
 prior={**v7['source_hashes'],**v6['source_hashes'],**v6['data_hashes'],'results/v6/headroom_results.json':v6['headroom_result_sha256'],**v5['source_hashes'],**v5['result_hashes']}
 preserved={p:sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in prior.items()}
 if not all(preserved.values()):raise AssertionError('historical source drift')
 sources=list((ROOT/'experiments').glob('v9*.py'))+list((ROOT/'tests').glob('test_v9*.py'))
 report=dict(status='verified',accepted_certificates_replayed=len(receipts),per_arm=dict(Counter(r['arm'] for r in receipts)),receipts=receipts,
             elapsed_seconds=perf_counter()-start,historical_hash_checks=preserved,v7_v8_archives_unchanged=True,v8_sources_unchanged=True,
             archive_sha256=sha256((ROOT/'results/v9/certificate_archive.json').read_bytes()).hexdigest(),
             source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
             acquisition='not_run',heldout_evaluation='not_run',compounding_capability='not_established',physical_validation='not_performed')
 (ROOT/'results/v9/verification_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print({k:report[k] for k in ('status','accepted_certificates_replayed','per_arm','elapsed_seconds')})
if __name__=='__main__':run()
