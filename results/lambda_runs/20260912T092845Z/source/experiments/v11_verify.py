"""Replay V11 saved certificates, compare exact outputs, preserve prior versions."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
from time import perf_counter
import json,re,platform
from .v7_verify import replay
from .v7_certificate import evaluate,clean
ROOT=Path(__file__).resolve().parents[1]


def run():
    start=perf_counter();receipts=[];signatures={}
    for domain,result_name,archive_name in [('matrix','preflight.json','certificate_archive.json'),
                                           ('guard','gate_cost.json','gate_certificate_archive.json')]:
        result=json.loads((ROOT/'results/v11'/result_name).read_text())
        archive=json.loads((ROOT/'results/v11'/archive_name).read_text())
        if len(archive)!=(2 if domain=='matrix' else 104):raise AssertionError('archive count')
        seen=set()
        for record in archive:
            arm=record['tag']['arm'];key=(record['n'],record['tag']['family'],record['gamma'],record['time'])
            if (arm,key) in seen:raise AssertionError('duplicate archive member')
            seen.add((arm,key));ans=replay(record)
            if ans['status']!='certified':raise AssertionError(ans)
            cs=tuple({p:F(v) for p,v in c.items()} for c in record['pieces'][0]['coefficients'])
            out=clean(evaluate(cs,F(record['time'])),record['n'],512)
            if out!={p:F(v) for p,v in record['output'].items()}:raise AssertionError('endpoint mismatch')
            signature=(cs,record['witness']['witnesses'],ans['bound'],out)
            if key in signatures and signature!=signatures[key]:raise AssertionError('exact cross-arm mismatch')
            signatures[key]=signature
            receipts.append(dict(domain=domain,arm=arm,n=key[0],family=key[1],gamma=key[2],time=key[3],bound=ans['bound']))
        for p,h in result['source_sha256'].items():
            if sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise AssertionError('V11 source drift '+p)
    for version in ('v8','v9','v10'):
        old=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
        for p,h in old['source_sha256'].items():
            if sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise AssertionError('prior source drift '+p)
    for version in ('v7','v8','v9','v10'):
        old=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
        if sha256((ROOT/f'results/{version}/certificate_archive.json').read_bytes()).hexdigest()!=old['archive_sha256']:raise AssertionError('prior archive drift')
    v7=json.loads((ROOT/'results/v7/verification_receipt.json').read_text());v6=json.loads((ROOT/'results/v6/verification_receipt.json').read_text());v5=json.loads((ROOT/'results/v5/verification_receipt.json').read_text())
    old={**v7['source_hashes'],**v6['source_hashes'],**v6['data_hashes'],'results/v6/headroom_results.json':v6['headroom_result_sha256'],**v5['source_hashes'],**v5['result_hashes']}
    preserved={p:sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old.items()}
    if not all(preserved.values()):raise AssertionError('historical drift')
    log=(ROOT/'results/v11/tests.log').read_text();count=re.findall(r'Ran (\d+) tests',log)
    if not count or not log.rstrip().endswith('OK'):raise AssertionError('test suite failed')
    sources=list((ROOT/'experiments').glob('v11*.py'))+list((ROOT/'tests').glob('test_v11*.py'))
    report=dict(status='verified',accepted_certificates_replayed=len(receipts),per_arm=dict(Counter(r['arm'] for r in receipts)),
                exact_cross_arm_calculations=len(signatures),receipts=receipts,elapsed_seconds=perf_counter()-start,
                historical_hash_checks=preserved,v7_v8_v9_v10_archives_unchanged=True,v8_v9_v10_sources_unchanged=True,
                archive_sha256={name:sha256((ROOT/'results/v11'/name).read_bytes()).hexdigest() for name in ('certificate_archive.json','gate_certificate_archive.json')},
                source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
                test_suite_count=int(count[-1]),test_suite_passed=True,test_log_sha256=sha256(log.encode()).hexdigest(),
                runtime=platform.python_version(),acquisition='not_run',heldout_evaluation='not_run',
                compounding_capability='not_established',physical_validation='not_performed')
    (ROOT/'results/v11/verification_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print({k:report[k] for k in ('status','accepted_certificates_replayed','per_arm','test_suite_count','elapsed_seconds')},flush=True)


if __name__=='__main__':run()
