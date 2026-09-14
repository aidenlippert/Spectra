"""Collect new exact pair-family caps and independent consistency probes."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/marginal_graded_hubbard8/pair_transfer'


def read(p):return json.loads(p.read_text())
def rel(p):return str(p.relative_to(ROOT))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n')


def main():
    log=(BASE/'family_full_validation.log').read_text()
    passed=re.search(r'(\d+) passed',log);subtests=re.search(r'(\d+) subtests passed',log)
    if not passed or ' failed' in log:raise ValueError('Completed passing full suite required')
    suites=list(ET.parse(BASE/'family_full_validation.xml').getroot().iter('testsuite'))
    if not suites or any(int(s.attrib.get(k,0)) for s in suites for k in ('failures','errors')):
        raise ValueError('Passing JUnit record required')
    if '30 passed' not in (BASE/'family_focused_validation.log').read_text():raise ValueError('Focused family checks missing')
    if '2 passed' not in (BASE/'two_spectator_validation.log').read_text():raise ValueError('Independent bit-action tests missing')
    names=['range_two_replay','range_two_family_limit_replay','pair_transfer_overlap',
           'one_spectator_overlap','spin_overlap','coherent_overlap','charge_markov_extension','two_spectator_overlap']
    paths=[BASE/c/'family'/(n+'.json') for c in ('W_zero','W_plus_1') for n in names]
    source_files=set();hash_entries=0
    for p in paths:
        d=read(p)
        if d.get('accepted') is not True:raise ValueError('Unaccepted current receipt: '+rel(p))
        for name,h in d['source_sha256'].items():
            if sha(ROOT/name)!=h:raise ValueError('Stale current source: '+name)
            hash_entries+=1;source_files.add(ROOT/name)
    cases={}
    for case in ('W_zero','W_plus_1'):
        p=BASE/case/'family';family=read(p/'range_two_family_limit_replay.json')
        energy=read(p/'range_two_replay.json');pair=read(p/'pair_transfer_overlap.json')
        two=read(p/'two_spectator_overlap.json')
        if not family.get('pair_transfer') or any(F(item['exact_moment']) for item in pair['separators']):
            raise ValueError('Exact pair closure missing')
        if two['new_independent_directions']!=30 or two['violated_directions']!=30:
            raise ValueError('Expected independently checked two-spectator obstruction missing')
        gap=F(family['periodic_family_upper'])-F(family['accepted_periodic_lower'])
        if gap!=F(family['family_gap']) or gap<0:raise ValueError('Exact family interval inconsistent')
        if F(energy['lower_replay']['periodic_lower_density'])!=F(family['accepted_periodic_lower']):
            raise ValueError('Energy/family lower mismatch')
        largest=max(two['moments'],key=lambda d:abs(F(d['exact_moment'])))
        cases[case]={'directory':rel(p),'periodic_lower':family['accepted_periodic_lower'],
                     'family_ceiling':family['periodic_family_upper'],'family_gap':str(gap),'family_gap_float':float(gap),
                     'mixture_sources':family['family_replay']['nearest_constraint_replay']['mixture_sources'],
                     'open_lower_per_site':energy['lower_per_site'],'physical_upper_per_site':energy['upper_per_site'],
                     'two_spectator_violations':30,'new_independent_directions':30,'largest_two_spectator_violation':largest}
    # Historical source changes remain visible; do not overwrite old receipts.
    historical=read(BASE/'summary.json')['current_receipts']+read(ROOT/'results/marginal_graded_hubbard8/spectator_hopping/thermal_limit_summary.json')['all_receipts_audited']
    snapshots={};snapshot_manifests=[]
    for directory in (BASE/'source_before',BASE/'family_source_before'):
        mp=directory/'original_sha256.json';snapshot_manifests.append(mp)
        for source,h in read(mp).items():
            p=directory/Path(source).name
            if sha(p)!=h:raise ValueError('Invalid preserved source snapshot')
            snapshots[(source,h)]=rel(p)
    unchanged=archived=0
    for name in historical:
        d=read(ROOT/name)
        if d.get('accepted') is not True:raise ValueError('Unaccepted historical receipt')
        for source,h in d['source_sha256'].items():
            if sha(ROOT/source)==h:unchanged+=1
            elif (source,h) in snapshots:archived+=1
            else:raise ValueError('Unaccounted historical source change: '+source)
    report=ROOT/'research/marginal_pair_family_limit.md'
    files=set(paths)|source_files|set(snapshot_manifests)|{report,Path(__file__).resolve()}
    for case in cases:files.update(p for p in (BASE/case/'family').rglob('*') if p.is_file())
    for directory in (BASE/'source_before',BASE/'family_source_before'):
        files.update(p for p in directory.rglob('*') if p.is_file())
    files.update(ROOT/name for name in historical)
    files.update([BASE/'summary.json',ROOT/'results/marginal_graded_hubbard8/spectator_hopping/thermal_limit_summary.json'])
    for case in cases:
        for path in (BASE/case/'family').rglob('*.json'):
            for name,h in read(path).get('source_sha256',{}).items():
                if sha(ROOT/name)!=h:raise ValueError('Stale construction input: '+name)
                files.add(ROOT/name)
    for name in ('pair_transfer_family_numeric.py','pair_transfer_thermal_atoms.py','two_spectator_overlap_probe.py',
                 'signed_charge_numeric.py','joint_profile_numeric.py'):
        files.add(ROOT/'results/marginal_graded_hubbard8/discovery'/name)
    files.update(ROOT/'tests'/name for name in ('test_marginal_pair_family.py','test_marginal_two_spectator_probe.py'))
    files.update(BASE/name for name in ('family_focused_validation.log','two_spectator_validation.log',
                                      'family_full_validation.log','family_full_validation.xml'))
    write(BASE/'family_provenance.json',{'scope':'Construction and current proof provenance. Numerical atoms and LP diagnostics remain nonaccepting. Historical source snapshots are retained separately from current replay claims.',
                                       'sha256':{rel(p):sha(p) for p in sorted(files)}})
    summary={'accepted':True,'cases':cases,'current_receipts':list(map(rel,paths)),
             'current_receipts_audited':len(paths),'current_source_hash_entries_checked':hash_entries,
             'historical_audit':{'receipts':len(historical),'unchanged_hash_entries':unchanged,
                                 'preserved_snapshot_entries':archived,'receipt_paths':historical,
                                 'scope':'Historical evidence checked against preserved original sources where current code changed. Not a fresh replay under the new verifier.'},
             'full_validation':{'command':'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --junitxml=results/marginal_graded_hubbard8/pair_transfer/family_full_validation.xml',
                                'tests_passed':int(passed[1]),'subtests_passed':int(subtests[1]) if subtests else 0,
                                'junit_suites':[s.attrib for s in suites],
                                'scope':'Full suite collected after FAMILY v10 changes. Two later two-spectator tests ran separately.'},
             'focused_family_tests_passed':30,'additional_two_spectator_tests_passed':2,
             'energy_verifier_changed':False,'family_version':10,'source_cap':89,
             'report':rel(report),'provenance':rel(BASE/'family_provenance.json'),
             'remaining':['No exact optimum attainment claim; certified finite family gaps are recorded above.',
                          'All30 independent two-spectator moments violate stationarity in both new mixtures; energy integration remains open.',
                          'No new transfer this turn. General quantum representability, generic molecular/long-range/higher-dimensional transfer and requested-accuracy scalability remain unproved. Goal active.']}
    write(BASE/'family_summary.json',summary)
    central=ROOT/'results/marginal_final_validation.json';d=read(central)
    d['pair_transfer_family_limit']=summary
    d['latest_validation_scope']='Exact81-source FAMILY v10 caps bracket W0/W1 within1.59e-6/2.02e-6 per site. Independent pair, spectator, spin and hopping closures pass. All30 independent two-spectator conditions fail in both mixtures. Current receipts replayed and full/focused regressions pass; historical source snapshots distinguished. No new energy improvement or transfer; general representability remains open.'
    write(central,d)
    print(json.dumps({'accepted':True,'current_receipts':len(paths),'current_hash_entries':hash_entries,
                      'historical_receipts':len(historical),'historical_snapshot_entries':archived,
                      'full_tests':int(passed[1]),'subtests':summary['full_validation']['subtests_passed'],
                      'provenance_files':len(files)}))


if __name__=='__main__':main()
