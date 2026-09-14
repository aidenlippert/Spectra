"""Collect verified thermal-limit receipts without treating discovery as proof."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spectator_hopping'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    previous = read(BASE / 'basis_completion_summary.json')
    names = ['range_two_replay', 'range_two_family_limit_replay', 'interval_refinement',
             'pair_transfer_overlap', 'one_spectator_overlap', 'spin_overlap',
             'coherent_overlap', 'charge_markov_extension']
    new_paths = [BASE / c / 'thermal/final' / (n + '.json')
                 for c in ('W_zero', 'W_plus_1') for n in names]
    all_paths = list(dict.fromkeys(previous['all_receipts_audited'] + list(map(relative, new_paths))))
    checked = 0
    for name in all_paths:
        receipt = read(ROOT / name)
        if receipt.get('accepted') is not True:
            raise ValueError('Receipt not accepted: ' + name)
        for source, expected in receipt['source_sha256'].items():
            if digest(ROOT / source) != expected:
                raise ValueError('Source changed: ' + source)
            checked += 1
    cases = {}
    for case in ('W_zero', 'W_plus_1'):
        folder = BASE / case / 'thermal/final'
        energy, family, refinement, pair = [read(folder / (n + '.json'))
            for n in ('range_two_replay', 'range_two_family_limit_replay',
                      'interval_refinement', 'pair_transfer_overlap')]
        if F(energy['lower_replay']['periodic_lower_density']) != F(family['accepted_periodic_lower']):
            raise ValueError('Energy/family lower mismatch')
        if not all(x['violated'] and F(x['exact_moment']) != 0 for x in pair['separators']):
            raise ValueError('Expected pair-transfer obstructions absent')
        cases[case] = {
            'directory': relative(folder), 'periodic_lower': family['accepted_periodic_lower'],
            'family_ceiling': family['periodic_family_upper'], 'family_gap': family['family_gap'],
            'family_gap_float': float(F(family['family_gap'])),
            'gap_reduction_factor': refinement['gap_reduction_factor'],
            'open_lower_per_site': energy['lower_per_site'],
            'physical_upper_per_site': energy['upper_per_site'],
            'pair_transfer_moments': pair['separators'],
            'lower_origin': 'thermal/scaled' if case == 'W_zero' else 'thermal/fullmetric',
            'family_origin': 'thermal/affine_completion' if case == 'W_zero' else 'thermal/final',
            'mixture_sources': len(read(folder / 'diagonal_family_limit_proposal.json')['mixture']),
        }
    log = BASE / 'thermal_focused_validation.log'
    if '25 passed in 0.07s' not in log.read_text():
        raise ValueError('Expected completed focused-test record absent')
    tests = [ROOT / 'tests/test_marginal_interval_refinement.py',
             ROOT / 'tests/test_marginal_fraction_free_completion.py']
    report = ROOT / 'research/marginal_spectator_thermal_limit.md'
    provenance = BASE / 'thermal_limit_provenance.json'
    files = set(new_paths + tests + [log, report, BASE / 'basis_completion_summary.json',
                                   BASE / 'basis_completion_provenance.json'])
    # Conservatively include all local discovery modules and physical modules,
    # including imports not listed by older proposal-only diagnostics.
    files.update((ROOT / 'results/marginal_graded_hubbard8/discovery').glob('*.py'))
    files.update((ROOT / 'experiments').glob('marginal*.py'))
    for case in cases:
        files.update(p for p in (BASE / case / 'thermal').rglob('*')
                     if p.is_file() and '__pycache__' not in p.parts)
    write(provenance, {'scope': 'Construction provenance, not physical acceptance. Includes failed proposals, numerical histories and all local discovery/physical modules. Final production receipts alone establish physical and family acceptance.',
                       'sha256': {relative(p): digest(p) for p in sorted(files)}})
    summary = {
        'accepted': True, 'cases': cases, 'new_receipts': list(map(relative, new_paths)),
        'all_receipts_audited': all_paths, 'receipts_audited': len(all_paths),
        'source_hash_entries_checked': checked,
        'focused_validation': {'command': 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q tests/test_marginal_interval_refinement.py tests/test_marginal_fraction_free_completion.py',
                               'tests_passed': 25, 'elapsed_seconds': 0.07,
                               'source_sha256': {relative(p): digest(p) for p in tests + [log]}},
        'production_verifier_changed': False, 'full_suite_rerun': False,
        'report': relative(report), 'provenance': relative(provenance),
        'scope': 'Exact brackets of width below 2e-7 per site for the optimum in the fixed W0/W1 spectator family. No exact attainment or physical ground-energy upper interpretation for family ceilings.',
        'remaining': ['All four pair-transfer consistency violations persist in both new mixtures; energy integration remains undone.',
                      'W=+/-0.1 retain older polished checkpoints; no new generic molecular or other-model transfer.',
                      'General quantum representability and requested-accuracy scalability remain unproved. Goal active.'],
    }
    write(BASE / 'thermal_limit_summary.json', summary)
    central = ROOT / 'results/marginal_final_validation.json'
    data = read(central)
    data['spectator_thermal_limit'] = summary
    data['latest_validation_scope'] = 'Thermal spectator limit: W0/W1 fixed-family gaps <2e-7/site; 16 new independently accepted receipts; 97 receipts audited; 25 focused tests. Production verifiers unchanged; full suite not rerun; pair-transfer energy integration and general representability remain open.'
    write(central, data)
    print(json.dumps({'accepted': True, 'receipts': len(all_paths), 'source_hash_entries': checked,
                      'provenance_files': len(files), 'gaps': {c: v['family_gap_float'] for c,v in cases.items()}}))


if __name__ == '__main__':
    main()
