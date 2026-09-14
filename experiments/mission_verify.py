"""Replay saved mission/V12 artifacts and check frozen historical hashes."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json,re
from .mission_policy_cover import Sample,certify
from .mission_examples import reaction,thermal,pack
from .v7_certificate import Generator,Piece,check_certificate,evaluate,clean
from .v8_fractional_cover import check_evolution_cover
ROOT=Path(__file__).resolve().parents[1]


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    examples=json.loads((ROOT/'results/mission/mathematical_examples.json').read_text())
    replays=[]
    for row in examples['rows']:
        result=row['result']
        samples=[Sample(tuple(F(v) for v in s['point']),tuple(tuple(F(v) for v in p) for p in s['intervals'])) for s in result['samples']]
        answer=certify(result['paths'],samples,[F(v) for v in row['lipschitz']],row['dimension'],F(row['epsilon']))
        if pack(answer)!=result['receipt']:raise AssertionError('saved policy certificate drift')
        physical=reaction if row['name']=='two_state_reaction' else thermal
        for sample in samples:
            values,_=physical(sample.point)
            if any(not lo<=v<=hi for v,(lo,hi) in zip(values,sample.intervals)):
                raise AssertionError('saved mathematical observation drift')
        replays.append(dict(name=row['name'],status=answer['status'],observations=len(samples),gap=str(answer['gap'])))
    v12=json.loads((ROOT/'results/v12/preflight.json').read_text())
    archive=json.loads((ROOT/'results/v12/certificate_archive.json').read_text());v12receipts=[]
    for record in archive:
        g=Generator({p:F(v) for p,v in record['h'].items()},F(record['gamma']),record['n'],512)
        initial={p:F(v) for p,v in record['initial'].items()}
        pieces=[Piece(F(p['duration']),tuple({k:F(v) for k,v in c.items()} for c in p['coefficients'])) for p in record['pieces']]
        checker=check_evolution_cover if record['witness']['schema']=='v8-evolution-cover-1' else check_certificate
        answer=checker(g,initial,pieces,record['witness'],F(record['tolerance']),expected_time=F(record['time']))
        if answer['status']!='certified':raise AssertionError(answer)
        out=clean(evaluate(pieces[0].coefficients,F(record['time'])),record['n'],512)
        if {k:str(v) for k,v in out.items()}!=record['output']:raise AssertionError('V12 endpoint drift')
        v12receipts.append(dict(arm=record['tag']['arm'],bound=answer['bound'],degree=len(pieces[0].coefficients)-1))
    if len(archive)!=3:raise AssertionError('V12 archive count')
    for p,h in {**examples['source_sha256'],**v12['source_sha256']}.items():
        if digest(p)!=h:raise AssertionError('new artifact source drift: '+p)
    for version in ('v8','v9','v10','v11'):
        old=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
        for p,h in old['source_sha256'].items():
            if digest(p)!=h:raise AssertionError('historical source drift: '+p)
    for version in ('v7','v8','v9','v10','v11'):
        old=json.loads((ROOT/f'results/{version}/verification_receipt.json').read_text())
        hashes=old['archive_sha256']
        if isinstance(hashes,str):hashes={'certificate_archive.json':hashes}
        for name,h in hashes.items():
            if digest(f'results/{version}/{name}')!=h:raise AssertionError('historical archive drift')
    v7=json.loads((ROOT/'results/v7/verification_receipt.json').read_text())
    v6=json.loads((ROOT/'results/v6/verification_receipt.json').read_text())
    v5=json.loads((ROOT/'results/v5/verification_receipt.json').read_text())
    historical={**v7['source_hashes'],**v6['source_hashes'],**v6['data_hashes'],
        'results/v6/headroom_results.json':v6['headroom_result_sha256'],**v5['source_hashes'],**v5['result_hashes']}
    checks={p:digest(p)==h for p,h in historical.items()}
    if not all(checks.values()):raise AssertionError('prior source/data/result drift')
    log=(ROOT/'results/mission/tests.log').read_text();counts=re.findall(r'Ran (\d+) tests',log)
    if not counts or not log.rstrip().endswith('OK'):raise AssertionError('full test suite did not pass')
    sources=list((ROOT/'experiments').glob('mission*.py'))+list((ROOT/'tests').glob('test_mission*.py'))
    report=dict(status='verified',mathematical_policy_replays=replays,v12_replays=v12receipts,
        v12_headroom_gate=v12['headroom_gate'],historical_hash_checks=checks,
        v7_through_v11_archives_unchanged=True,v8_through_v11_sources_unchanged=True,
        tests_passed=int(counts[-1]),test_log_sha256=sha256(log.encode()).hexdigest(),
        source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
        physical_validation='not_performed',autonomous_acquisition='not_performed',
        general_matter_mastery='not_established')
    (ROOT/'results/mission/verification_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    v12report=dict(status='verified',certificates_replayed=3,receipts=v12receipts,
        headroom_gate=v12['headroom_gate'],archive_sha256=digest('results/v12/certificate_archive.json'),
        historical_hash_checks=checks,v7_through_v11_archives_unchanged=True,
        source_sha256={p:digest(p) for p in ['experiments/v12_preflight.py','experiments/v12_small_cover.py',
            'experiments/v12_small_cover_taylor.py','tests/test_v12_small_cover.py']},
        source_scope='complete-cost preflight implementation and focused tests; diagnostic-only helper scripts are not dependencies',
        test_suite_count=int(counts[-1]),test_suite_passed=True,physical_validation='not_performed',acquisition='not_run')
    (ROOT/'results/v12/verification_receipt.json').write_text(json.dumps(v12report,indent=2)+'\n')
    print(json.dumps(dict(status='verified',mathematical_policy_replays=replays,v12_replays=v12receipts,tests_passed=int(counts[-1]),historical_checks=len(checks)),indent=2))


if __name__=='__main__':run()
