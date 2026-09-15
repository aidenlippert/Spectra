"""Replay all saved V10 calculations and verify historical preservation."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from time import perf_counter
from collections import Counter
import json,re,platform
from .v7_verify import replay
from .v7_certificate import Piece,evaluate,clean
from .v10_headroom import replay_reference,model
ROOT=Path(__file__).resolve().parents[1]


def run():
    start=perf_counter();saved=json.loads((ROOT/'results/v10/headroom.json').read_text())
    archive=json.loads((ROOT/'results/v10/certificate_archive.json').read_text())
    expected={(r['arm'],r['n'],r['family'],r['gamma'],r['time']):r for r in saved['rows'] if r['status']=='certified'}
    if len(expected)!=len(archive):raise AssertionError('archive count mismatch')
    receipts=[]
    for record in archive:
        if record.get('schema')=='v10-reference-1':
            result=replay_reference(record);arm='reference'
            families=[family for family in ('xxz','mixed') if {p:str(c) for p,c in model(record['n'],family)[0].items()}==record['h']]
            if len(families)!=1:raise AssertionError('ambiguous physical workload')
            family=families[0]
        else:
            arm,family=record['tag']['arm'],record['tag']['family'];result=replay(record)
            cs=tuple({p:F(c) for p,c in d.items()} for d in record['pieces'][0]['coefficients'])
            output=clean(evaluate(cs,F(record['time'])),record['n'],512)
            if output!={p:F(c) for p,c in record['output'].items()}:raise AssertionError('polynomial output mismatch')
        key=(arm,record['n'],family,record['gamma'],record['time'])
        if key not in expected:raise AssertionError('duplicate/unknown archive member')
        row=expected.pop(key)
        if result['status']!='certified' or result['bound']!=row['bound']:raise AssertionError('bound mismatch')
        receipts.append(dict(arm=arm,n=record['n'],family=family,gamma=record['gamma'],time=record['time'],bound=result['bound']))
    if expected:raise AssertionError('missing archive members')
    for p,h in saved['source_sha256'].items():
        if sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise AssertionError('V10 experiment source drift '+p)
    for version in ('v8','v9'):
        prior=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
        for p,h in prior['source_sha256'].items():
            if sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise AssertionError('prior source drift '+p)
    for version in ('v7','v8','v9'):
        prior=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
        if sha256((ROOT/f'results/{version}/certificate_archive.json').read_bytes()).hexdigest()!=prior['archive_sha256']:raise AssertionError('prior archive drift')
    v7=json.loads((ROOT/'results/v7/verification_receipt.json').read_text());v6=json.loads((ROOT/'results/v6/verification_receipt.json').read_text());v5=json.loads((ROOT/'results/v5/verification_receipt.json').read_text())
    old={**v7['source_hashes'],**v6['source_hashes'],**v6['data_hashes'],
         'results/v6/headroom_results.json':v6['headroom_result_sha256'],**v5['source_hashes'],**v5['result_hashes']}
    preserved={p:sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old.items()}
    if not all(preserved.values()):raise AssertionError('historical drift')
    log=(ROOT/'results/v10/tests.log').read_text();counts=re.findall(r'Ran (\d+) tests',log)
    if not counts or not log.rstrip().endswith('OK'):raise AssertionError('test suite did not pass')
    sources=list((ROOT/'experiments').glob('v10*.py'))+list((ROOT/'tests').glob('test_v10*.py'))
    report=dict(status='verified',accepted_certificates_replayed=len(receipts),per_arm=dict(Counter(r['arm'] for r in receipts)),
                receipts=receipts,elapsed_seconds=perf_counter()-start,historical_hash_checks=preserved,
                v7_v8_v9_archives_unchanged=True,v8_v9_sources_unchanged=True,
                archive_sha256=sha256((ROOT/'results/v10/certificate_archive.json').read_bytes()).hexdigest(),
                archive_bytes=(ROOT/'results/v10/certificate_archive.json').stat().st_size,
                source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
                test_suite_count=int(counts[-1]),test_suite_passed=True,test_log_sha256=sha256(log.encode()).hexdigest(),
                runtime=platform.python_version(),acquisition='not_run',heldout_evaluation='not_run',
                compounding_capability='not_established',physical_validation='not_performed')
    (ROOT/'results/v10/verification_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print({k:report[k] for k in ('status','accepted_certificates_replayed','per_arm','test_suite_count','elapsed_seconds')},flush=True)


if __name__=='__main__':run()
