"""Collect exact full-overlap obstructions and preserve prior energy evidence."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'results/marginal_graded_hubbard8/full_overlap'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    prior_root = BASE.parent/'three_spectator'
    prior = read(prior_root/'summary.json')
    prior_manifest = read(prior_root/'provenance.json')['sha256']
    for name, expected in prior_manifest.items():
        if sha(ROOT/name) != expected:
            raise ValueError('Prior provenance changed: '+name)
    suites = list(ET.parse(BASE/'validation.xml').getroot().iter('testsuite'))
    if (len(suites) != 1 or int(suites[0].attrib['tests']) != 21
            or any(int(suites[0].attrib.get(k, 0)) for k in ('errors', 'failures', 'skipped'))
            or '21 passed' not in (BASE/'validation.log').read_text()):
        raise ValueError('Twenty-one passing focused tests required')
    receipts = [BASE/case/name for case in ('W_zero', 'W_plus_1')
                for name in ('full_overlap_replay.json', 'telescope_replay.json')]
    receipts.append(BASE/'support_replay.json')
    files = {Path(__file__).resolve(), prior_root/'summary.json', prior_root/'provenance.json',
             ROOT/'research/marginal_full_overlap_obstruction.md'}
    entries = 0
    for path in receipts:
        receipt = read(path)
        if receipt.get('accepted') is not True:
            raise ValueError('Unaccepted new proof receipt')
        for name, expected in receipt['source_sha256'].items():
            source = ROOT/name
            if sha(source) != expected:
                raise ValueError('New receipt source changed: '+name)
            files.add(source); entries += 1
    cases = {}
    for case in ('W_zero', 'W_plus_1'):
        overlap = read(BASE/case/'full_overlap_replay.json')
        telescope = read(BASE/case/'telescope_replay.json')
        w = overlap['witness']; delta = F(w['left_probability'])-F(w['right_probability'])
        if (not delta or delta != F(w['exact_difference'])
                or delta != F(telescope['exact_direct_mixture_expectation'])
                or not overlap['stationary_extension_refuted']
                or not 0 <= F(w['left_probability']) <= 1
                or not 0 <= F(w['right_probability']) <= 1
                or telescope['energy_certificate'] is not False):
            raise ValueError('Overlap, positive-projector and direct-moment receipts disagree')
        cases[case] = {
            'full_overlap_receipt': str((BASE/case/'full_overlap_replay.json').relative_to(ROOT)),
            'telescope_receipt': str((BASE/case/'telescope_replay.json').relative_to(ROOT)),
            'stationary_extension_of_this_mixture_refuted': True,
            'five_site_vector': w['five_site_vector'], 'exact_mismatch': str(delta),
            'mismatch_float': float(delta), 'difference_upper_nonzeros': overlap['difference_upper_nonzeros'],
            'difference_diagonal_nonzeros': overlap['difference_diagonal_nonzeros'],
            'common_denominator_bits': overlap['common_denominator_bits'],
            'five_site_telescope_nonzeros': telescope['five_site_nonzeros'],
            'six_site_telescope_nonzeros': telescope['six_site_nonzeros'],
            'operator_norm_bound': 1,
            'unchanged_periodic_lower': prior['cases'][case]['periodic_lower'],
            'unchanged_family_gap': prior['cases'][case]['family_gap']}
    support = read(BASE/'support_replay.json')
    if (support['prior_telescope_count'] != 71
            or support['independent_new_directions_modulo_prior_family'] != 2
            or not F(support['new_direction_minor_determinant'])):
        raise ValueError('Independent new directions are not proved')
    files.update(p for p in BASE.rglob('*') if p.is_file() and p.name not in ('summary.json', 'provenance.json'))
    files.update((BASE.parent/'discovery').glob('full_overlap*.py'))
    files.update((ROOT/'tests').glob('test_marginal_full_overlap*.py'))
    result = {
        'accepted': True, 'new_energy_certificate': False, 'goal_status': 'active',
        'cases': cases, 'independent_new_directions': 2,
        'current_receipts': [str(p.relative_to(ROOT)) for p in receipts],
        'current_receipts_audited': len(receipts), 'current_source_hash_entries_checked': entries,
        'prior_provenance_files_checked': len(prior_manifest),
        'prior_energy_summary_sha256': sha(prior_root/'summary.json'),
        'prior_energy_provenance_sha256': sha(prior_root/'provenance.json'),
        'focused_tests_passed': 21, 'focused_junit': suites[0].attrib,
        'production_sources_changed': False, 'full_suite_rerun': False,
        'prior_full_suite': {'tests_passed': 982, 'subtests_passed': 102,
                            'separate_later_fraction_free_tests': 3},
        'report': 'research/marginal_full_overlap_obstruction.md',
        'provenance': str((BASE/'provenance.json').relative_to(ROOT)),
        'remaining': ['Integrate both independent stationary constraints and require fresh energy PSD and family moment acceptance.',
                      'Current lower bounds and nonzero family gaps are unchanged; no numerical-limit resolution.',
                      'Only these particular symmetry-averaged mixture extensions are refuted; general representability remains open.',
                      'No new coupling, geometry or molecular transfer. Requested-accuracy scalability remains unproved.']}
    write(BASE/'provenance.json', {'sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(files)},
                                 'scope': 'Five exact overlap/support receipts, their sources, tests, report and preserved prior energy summary. No new energy bound.'})
    write(BASE/'summary.json', result)
    central_path = ROOT/'results/marginal_final_validation.json'; central = read(central_path)
    central['full_overlap_obstruction'] = result
    central['latest_validation_scope'] = ('Exact full five-site overlap mismatches refute stationary extensions of both selected three-spectator local mixtures. Two compact symmetry-compatible telescopes add independent directions beyond the prior71 offdiagonal constraints. Five receipts and prior479-file provenance audited;21 focused tests pass. Production and energy bounds unchanged; goal active.')
    write(central_path, central)
    print(json.dumps({'accepted': True, 'current_receipts': len(receipts), 'source_hash_entries': entries,
                      'prior_provenance_files': len(prior_manifest), 'new_provenance_files': len(files),
                      'focused_tests_passed': 21, 'new_energy_certificate': False}))


if __name__ == '__main__':
    main()
