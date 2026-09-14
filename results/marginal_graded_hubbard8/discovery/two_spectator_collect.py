"""Audit exact selected certificates separately from untrusted construction data."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/two_spectator'
OLD = BASE.parent / 'pair_transfer'


def read(p): return json.loads(p.read_text())
def rel(p): return str(p.relative_to(ROOT))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p, d): p.write_text(json.dumps(d, indent=2) + '\n')


def main():
    log = (BASE / 'full_validation.log').read_text()
    passed = re.search(r'(\d+) passed', log)
    subtests = re.search(r'(\d+) subtests passed', log)
    if not passed or ' failed' in log:
        raise ValueError('Completed passing full suite required')
    suites = list(ET.parse(BASE / 'full_validation.xml').getroot().iter('testsuite'))
    if not suites or any(int(s.attrib.get(k, 0)) for s in suites for k in ('failures', 'errors')):
        raise ValueError('Passing JUnit record required')
    if '44 passed' not in (BASE / 'focused_validation.log').read_text():
        raise ValueError('Focused production checks missing')
    if '1 passed' not in (BASE / 'matching_validation.log').read_text():
        raise ValueError('Matching refusal check missing')
    names = ['range_two_replay', 'range_two_family_limit_replay', 'strict_family_separation',
             'two_spectator_closure', 'pair_transfer_overlap', 'one_spectator_overlap',
             'spin_overlap', 'coherent_overlap', 'charge_markov_extension']
    paths = [BASE / c / 'final' / (n + '.json') for c in ('W_zero', 'W_plus_1') for n in names]
    paths += [BASE / c / 'previous_family' / (n + '.json')
              for c in ('W_zero', 'W_plus_1') for n in ('range_two_replay', 'range_two_family_limit_replay')]
    paths += [BASE / 'held_out' / d / 'range_two_replay.json' for d in ('', 'without_two', 'previous_recipe')]
    paths += [BASE / 'held_out/frozen_transfer_comparison.json']
    source_files = set()
    hash_entries = 0
    for p in paths:
        d = read(p)
        if d.get('accepted') is not True:
            raise ValueError('Unaccepted current receipt: ' + rel(p))
        for name, h in d['source_sha256'].items():
            if sha(ROOT / name) != h:
                raise ValueError('Stale current source: ' + name)
            hash_entries += 1
            source_files.add(ROOT / name)
    cases = {}
    for case in ('W_zero', 'W_plus_1'):
        p = BASE / case / 'final'
        family = read(p / 'range_two_family_limit_replay.json')
        energy = read(p / 'range_two_replay.json')
        separation = read(p / 'strict_family_separation.json')
        two = read(p / 'two_spectator_closure.json')['two_spectator_moments']
        if not family.get('two_spectator_hopping') or len(two) != 30 or any(map(F, two.values())):
            raise ValueError('Exact two-spectator closure missing')
        for n in ('pair_transfer_overlap', 'spin_overlap'):
            if any(F(item['exact_moment']) for item in read(p / (n + '.json'))['separators']):
                raise ValueError('Preceding closure failed: ' + n)
        if any(map(F, read(p / 'one_spectator_overlap.json')['spectator_moments'].values())):
            raise ValueError('One-spectator closure failed')
        if F(read(p / 'coherent_overlap.json')['exact_moment']):
            raise ValueError('Coherent closure failed')
        if read(p / 'charge_markov_extension.json')['classical_only'] is not True:
            raise ValueError('Classical extension scope missing')
        gap = F(family['periodic_family_upper']) - F(family['accepted_periodic_lower'])
        if gap != F(family['family_gap']) or gap < 0:
            raise ValueError('Exact family interval inconsistent')
        if F(energy['lower_replay']['periodic_lower_density']) != F(family['accepted_periodic_lower']):
            raise ValueError('Energy/family lower mismatch')
        if F(separation['exact_separation']) <= 0:
            raise ValueError('Strict separation missing')
        cases[case] = {'directory': rel(p), 'periodic_lower': family['accepted_periodic_lower'],
                       'family_ceiling': family['periodic_family_upper'], 'family_gap': str(gap),
                       'family_gap_float': float(gap),
                       'mixture_sources': family['family_replay']['nearest_constraint_replay']['mixture_sources'],
                       'open_lower_per_site': energy['lower_per_site'], 'physical_upper_per_site': energy['upper_per_site'],
                       'strict_excess_over_previous_family': separation['exact_separation'],
                       'strict_excess_float': separation['separation_float'], 'two_spectator_zero_moments': 30}
    prior = read(OLD / 'family_summary.json')
    historical = prior['current_receipts'] + prior['historical_audit']['receipt_paths']
    if len(historical) != len(set(historical)):
        raise ValueError('Duplicate historical receipts')
    snapshots = {}
    snapshot_dirs = (OLD / 'source_before', OLD / 'family_source_before', BASE / 'source_before')
    for directory in snapshot_dirs:
        for source, h in read(directory / 'original_sha256.json').items():
            p = directory / Path(source).name
            if sha(p) != h:
                raise ValueError('Invalid preserved source snapshot')
            snapshots[(source, h)] = rel(p)
    unchanged = archived = 0
    for name in historical:
        d = read(ROOT / name)
        if d.get('accepted') is not True:
            raise ValueError('Unaccepted historical receipt')
        for source, h in d['source_sha256'].items():
            if sha(ROOT / source) == h:
                unchanged += 1
            elif (source, h) in snapshots:
                archived += 1
            else:
                raise ValueError('Unaccounted historical source change: ' + source)
    report = ROOT / 'research/marginal_two_spectator_energy.md'
    files = set(paths) | source_files | {report, Path(__file__).resolve(), OLD / 'family_summary.json'}
    files.update(ROOT / name for name in historical)
    files.update(p for p in BASE.rglob('*') if p.is_file() and p not in (BASE / 'summary.json', BASE / 'provenance.json'))
    for directory in snapshot_dirs:
        files.update(p for p in directory.rglob('*') if p.is_file())
    construction_entries = 0
    # Hash validity is provenance only: proposals do not become accepted proofs.
    for path in BASE.rglob('*.json'):
        if path.name in ('summary.json', 'provenance.json'):
            continue
        d = read(path)
        if not isinstance(d, dict):
            continue
        for name, h in d.get('source_sha256', {}).items():
            if sha(ROOT / name) != h:
                raise ValueError('Stale construction input: ' + name)
            construction_entries += 1
            files.add(ROOT / name)
    discovery = BASE.parent / 'discovery'
    files.update(discovery.glob('two_spectator*.py'))
    for name in ('joint_congruence_replay.py', 'range_two_density_replay.py', 'range_two_family_limit_replay.py',
                 'signed_charge_numeric.py', 'joint_profile_numeric.py', 'pair_transfer_separation.py',
                 'pair_transfer_target_replay.py', 'pair_transfer_probe.py', 'spectator_overlap_replay.py',
                 'spin_overlap_probe.py', 'coherent_overlap_probe.py', 'charge_markov_extension.py'):
        files.add(discovery / name)
    files.update(ROOT / 'tests' / n for n in ('test_marginal_two_spectator_hopping.py', 'test_marginal_two_spectator_matching.py'))
    diagnostics = {}
    for path in BASE.rglob('diagonal_family_limit_diagnostic.json'):
        d = read(path)
        diagnostics[rel(path)] = {k: d[k] for k in ('accepted', 'proposal_written', 'states', 'seconds', 'mixture_sources', 'rational_family_upper_float')}
        diagnostics[rel(path)]['pricing_rounds'] = len(d['pricing'])
        diagnostics[rel(path)]['last_pricing'] = d['pricing'][-1:] or []
    write(BASE / 'provenance.json', {'scope': 'Current construction and proof provenance. Hashing numerical proposals does not confer acceptance. Historical receipts are checked against preserved original source snapshots.',
                                      'sha256': {rel(p): sha(p) for p in sorted(files)}})
    summary = {'accepted': True, 'cases': cases, 'current_receipts': list(map(rel, paths)),
               'current_receipts_audited': len(paths), 'current_source_hash_entries_checked': hash_entries,
               'all_construction_source_hash_entries_checked': construction_entries,
               'historical_audit': {'receipts': len(historical), 'unchanged_hash_entries': unchanged,
                                    'preserved_snapshot_entries': archived, 'receipt_paths': historical,
                                    'scope': 'Historical evidence audited against preserved sources; not fresh replay under current verifier.'},
               'frozen_transfer': read(BASE / 'held_out/frozen_transfer_comparison.json'),
               'full_validation': {'command': 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --junitxml=results/marginal_graded_hubbard8/two_spectator/full_validation.xml',
                                   'tests_passed': int(passed[1]), 'subtests_passed': int(subtests[1]),
                                   'junit_suites': [s.attrib for s in suites], 'elapsed_seconds': 849.76,
                                   'scope': 'Final full suite includes both new production and matching tests. No subsequent production edits.'},
               'focused_tests_passed': 44, 'matching_tests_passed': 1,
               'energy_version': 15, 'family_version': 11, 'source_cap': 119,
               'discovery_diagnostics': diagnostics,
               'report': rel(report), 'provenance': rel(BASE / 'provenance.json'),
               'remaining': ['W=1 family gap remains 0.005548129720482485; no tight numerical-limit claim.',
                             'W=0 finite family gap is 0.000019009152447824; no exact attainment claim.',
                             'New mixtures pass tested necessary conditions; full quantum extension remains unproved.',
                             'Finite coupling transfer only; general molecular, long-range, higher-dimensional transfer and requested-accuracy scalability remain open. Goal active.']}
    write(BASE / 'summary.json', summary)
    central = ROOT / 'results/marginal_final_validation.json'
    d = read(central)
    d['two_spectator_energy'] = summary
    d['latest_validation_scope'] = 'ENERGY v15 and FAMILY v11 integrate all30 two-spectator conditions. Exact W0/W1 lower bounds strictly exceed preceding pair-family ceilings. Selected family gaps 1.901e-5 and0.005549; W1 limit remains unresolved. All tested closures pass. Frozen coupling transfer and ablation accepted. Full945 tests and102 subtests pass. Current sources and historical snapshots audited separately; general representability and requested-accuracy scalability open.'
    write(central, d)
    print(json.dumps({'accepted': True, 'current_receipts': len(paths), 'current_hash_entries': hash_entries,
                      'historical_receipts': len(historical), 'historical_unchanged_entries': unchanged,
                      'historical_snapshot_entries': archived, 'construction_entries': construction_entries,
                      'full_tests': int(passed[1]), 'subtests': int(subtests[1]), 'provenance_files': len(files)}))


if __name__ == '__main__':
    main()
