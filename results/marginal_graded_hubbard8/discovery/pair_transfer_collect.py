"""Finalize current-source pair-transfer evidence after completed regression."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'results/marginal_graded_hubbard8/pair_transfer'


def read(p): return json.loads(p.read_text())
def rel(p): return str(p.relative_to(ROOT))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,data): p.write_text(json.dumps(data,indent=2)+'\n')


def main():
    full_log = BASE/'full_validation.log'; log = full_log.read_text()
    passed = re.search(r'(\d+) passed',log)
    subtests = re.search(r'(\d+) subtests passed',log)
    if not passed or ' failed' in log: raise ValueError('Completed passing full suite required')
    tree = ET.parse(BASE/'full_validation.xml')
    suites = list(tree.getroot().iter('testsuite'))
    if not suites or any(int(s.attrib.get(k,0)) for s in suites for k in ('failures','errors')):
        raise ValueError('Passing JUnit suites required')
    selection = {'W_zero':BASE/'W_zero','W_plus_1':BASE/'W_plus_1/polished'}
    paths = []
    for directory in list(selection.values())+[BASE/'W_plus_1']:
        paths.extend(directory/name for name in ['range_two_replay.json','strict_family_separation.json',
                     'previous_family/range_two_replay.json','previous_family/range_two_family_limit_replay.json'])
    paths += [BASE/name for name in ['held_out/range_two_replay.json',
               'held_out/without_pair/range_two_replay.json','held_out/frozen_transfer_comparison.json']]
    paths += [directory/'previous_family/pair_transfer_overlap.json' for directory in selection.values()]
    paths = list(dict.fromkeys(paths)); hash_count = 0; proof_sources = set()
    for path in paths:
        d = read(path)
        if d.get('accepted') is not True: raise ValueError('Unaccepted receipt '+rel(path))
        for name,expected in d['source_sha256'].items():
            if sha(ROOT/name) != expected: raise ValueError('Stale current proof: '+name)
            hash_count += 1; proof_sources.add(ROOT/name)
    cases = {}
    for case,directory in selection.items():
        e = read(directory/'range_two_replay.json'); s = read(directory/'strict_family_separation.json')
        if F(s['exact_separation']) <= 0: raise ValueError('Strict separation absent')
        cases[case] = {'directory':rel(directory),'periodic_lower':e['lower_replay']['periodic_lower_density'],
                       'open_lower_per_site':e['lower_per_site'],'physical_upper_per_site':e['upper_per_site'],
                       'old_family_ceiling':s['preceding_family_ceiling'],
                       'strict_separation':s['exact_separation'],'strict_separation_float':s['separation_float'],
                       'pair_coefficients':e['lower_replay']['pair_transfer']}
    focused_log = (BASE/'focused_final_validation.log').read_text()
    separation_log = (BASE/'separation_validation.log').read_text()
    if '67 passed' not in focused_log or '9 passed' not in separation_log:
        raise ValueError('Completed focused checks required')
    report = ROOT/'research/marginal_pair_transfer_energy.md'
    files = set(p for p in BASE.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                and p.name not in ('summary.json','provenance.json'))
    files.update(proof_sources)
    files.update((ROOT/'results/marginal_graded_hubbard8/discovery').glob('pair_transfer*.py'))
    files.update((ROOT/'tests').glob('test_marginal_pair*.py'))
    files.update([ROOT/'experiments/marginal_pair_transfer.py',report])
    # Include imported numerical chart sources omitted by early proposal manifests.
    for name in ('joint_profile_numeric.py','signed_charge_numeric.py'):
        files.add(ROOT/'results/marginal_graded_hubbard8/discovery'/name)
    write(BASE/'provenance.json',{'scope':'Artifact and source provenance. Includes rejected numerical improvements, pre-change source snapshots and current receipts. Numerical artifacts are not accepting proof gates.',
                                 'sha256':{rel(p):sha(p) for p in sorted(files)}})
    summary = {'accepted':True,'cases':cases,'held_out':read(BASE/'held_out/frozen_transfer_comparison.json'),
               'current_receipts':list(map(rel,paths)),'current_receipts_audited':len(paths),
               'current_source_hash_entries_checked':hash_count,
               'historical_source_audit':read(BASE/'historical_source_audit.json'),
               'production_changes':['ENERGY v14 pair-transfer extension with all preceding gates retained.',
                                     'Old energy versions and old family matching refuse unsupported pair coefficients.'],
               'full_validation':{'command':'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --junitxml=results/marginal_graded_hubbard8/pair_transfer/full_validation.xml',
                                  'tests_passed':int(passed[1]),'subtests_passed':int(subtests[1]) if subtests else 0,
                                  'junit_suites':[s.attrib for s in suites],
                                  'log':rel(full_log),'xml':rel(BASE/'full_validation.xml'),
                                  'scope':'Full suite collected after production changes. The subsequently added nine separation/driver tests ran separately.'},
               'focused_tests_passed':67,'additional_separation_tests_passed':9,
               'report':rel(report),'provenance':rel(BASE/'provenance.json'),
               'remaining':['No certified numerical ceiling yet for the enlarged pair-transfer family.',
                            'Transfer changes couplings on the same geometry; no generic molecular, long-range or higher-dimensional result.',
                            'General quantum representability and requested-accuracy scalability remain unproved. Goal active.']}
    write(BASE/'summary.json',summary)
    central = ROOT/'results/marginal_final_validation.json';data = read(central)
    data['pair_transfer_energy'] = summary
    data['latest_validation_scope'] = 'Exact pair-transfer lower certificates exceed the entire preceding fixed family for W0/W1. Frozen recipe transfers to U5,V1/4,W-1/5; exact pair ablation improves by0.0242373/site. Current proofs freshly replayed, full and focused regressions passed. Historical source snapshots are distinguished from current-code evidence; enlarged-family ceiling and general representability remain open.'
    write(central,data)
    print(json.dumps({'accepted':True,'current_receipts':len(paths),'current_hash_entries':hash_count,
                      'full_tests_passed':int(passed[1]),'subtests_passed':summary['full_validation']['subtests_passed'],
                      'provenance_files':len(files)}))


if __name__ == '__main__': main()
