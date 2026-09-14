"""Collect existing exact receipts; does not substitute for physical replay."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spectator_hopping'
CASES = ('W_zero', 'W_plus_1_10', 'W_minus_1_10', 'W_plus_1')
RECEIPTS = ('range_two_replay', 'range_two_family_limit_replay', 'family_escape',
            'spectator_ablation', 'pair_transfer_overlap', 'one_spectator_overlap',
            'spin_overlap', 'coherent_overlap', 'charge_markov_extension')


def read(path):
    return json.loads(path.read_text())


def relative(path):
    return str(path.relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    checks = 0
    receipt_paths = []
    cases = {}
    provenance = {}
    for case in CASES:
        parent = BASE / case
        directory = parent / 'polished'
        for path in [parent / 'basis_probe.json'] + [directory / (n + '.json') for n in RECEIPTS]:
            result = read(path)
            assert result['accepted'], path
            for name, expected in result['source_sha256'].items():
                assert digest(ROOT / name) == expected, (path, name)
                checks += 1
            receipt_paths.append(relative(path))
        energy = read(directory / 'range_two_replay.json')
        family = read(directory / 'range_two_family_limit_replay.json')
        escape = read(directory / 'family_escape.json')
        ablation = read(directory / 'spectator_ablation.json')
        pair = read(directory / 'pair_transfer_overlap.json')
        basis = read(parent / 'basis_probe.json')
        assert len(basis['separators']) == 14
        assert all(F(s['exact_moment']) for s in basis['separators'])
        assert basis['exact_rank']['rank'] == 17
        assert pair['exact_rank']['rank'] == 4
        assert all(s['violated'] and F(s['exact_moment']) for s in pair['separators'])
        assert all(F(v) == 0 for v in read(directory / 'one_spectator_overlap.json')['spectator_moments'].values())
        assert F(ablation['without_spectator_residual']) < 0 <= F(ablation['with_spectator_residual'])
        assert F(escape['strict_escape']) > 0
        mixture = read(directory / 'diagonal_family_limit_proposal.json')['mixture']
        assert len(mixture) == 85
        old = read(ROOT / 'results/marginal_graded_hubbard8/spin_telescope' / case / 'polished/range_two_replay.json')
        proposals = [read(parent / 'spectator_proposal.json'), read(directory / 'spectator_proposal.json')]
        assert all(not p['accepted'] and not p['optimizer_success'] for p in proposals)
        cases[case] = {
            'target': energy['target'], 'directory': relative(directory),
            'lower_per_site': energy['lower_per_site'], 'lower_float': float(F(energy['lower_per_site'])),
            'upper_per_site': energy['upper_per_site'], 'upper_float': float(F(energy['upper_per_site'])),
            'gain_over_spin_lower': str(F(energy['lower_per_site']) - F(old['lower_per_site'])),
            'strict_escape_over_spin_family': escape['strict_escape'],
            'strict_escape_float': escape['strict_escape_float'],
            'accepted_periodic_lower': family['accepted_periodic_lower'],
            'periodic_family_ceiling': family['periodic_family_upper'],
            'periodic_family_gap': family['family_gap'], 'family_gap_float': family['family_gap_float'],
            'family_sources': len(mixture),
            'primal_passes': [{k: p[k] for k in ('optimizer_success', 'optimizer_message', 'matrix_evaluations', 'seconds')} for p in proposals],
            'ablation_without': ablation['without_spectator_residual'],
            'ablation_with': ablation['with_spectator_residual'],
            'pair_transfer_separators': pair['separators'],
        }
        rerun = case in ('W_minus_1_10', 'W_plus_1')
        seed = directory if rerun else parent
        driver = ROOT / 'results/marginal_graded_hubbard8/discovery/spectator_family_numeric.py' if rerun else BASE / 'pre_spectator_sources/spectator_family_before_multieigen.py'
        diagnostic = seed / 'diagonal_family_limit_diagnostic.json'
        assert read(diagnostic)['proposal_written']
        assert digest(seed / 'diagonal_family_limit_proposal.json') == digest(directory / 'diagonal_family_limit_proposal.json')
        files = [driver, seed / 'profile_joint_r1_2_certificate.json', diagnostic,
                 directory / 'diagonal_family_limit_proposal.json']
        provenance[case] = {
            'seed_directory': relative(seed), 'preserved_driver': relative(driver),
            'pricing_eigenvectors': 2 if rerun else 1, 'pricing_rounds': 40 if rerun else 60,
            'determinant_anchors': True, 'scale': read(diagnostic)['scale'], 'pricing_radius': 0,
            'diagnostic': relative(diagnostic),
            'sha256': {relative(p): digest(p) for p in files},
            'scope': 'Historical numerical proposal provenance. Final acceptance uses current production replay. Seed and round count also changed for the two-eigenvector reruns; this is not a controlled performance comparison.',
        }
    xml = ET.parse(BASE / 'full_validation.xml')
    suite = xml.find('testsuite')
    assert suite is not None
    assert all(suite.attrib[k] == '0' for k in ('errors', 'failures', 'skipped'))
    assert len(xml.findall('.//testcase')) == 847 and int(suite.attrib['tests']) == 949
    validation = {
        'command': 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q --junitxml=results/marginal_graded_hubbard8/spectator_hopping/full_validation.xml',
        'tests_passed': 847, 'subtests_passed': 102, 'elapsed_seconds': 661.31,
        'junit_tests_including_subtests': 949, 'junit_testcase_elements': 847,
        'errors': 0, 'failures': 0, 'skipped': 0,
        'warning': 'Existing test_v4_calibration.py::test_world returns a tuple.',
        'log': relative(BASE / 'full_validation.log'), 'junit': relative(BASE / 'full_validation.xml'),
    }
    limitations = [
        'All eight primal optimization passes reached the iteration limit; the current 85-variable family optimum remains unresolved. The exact family gaps are loose.',
        'Strict escape covers the same fixed projector sources, ratio, ceilings and nine-shape sparse span, with prior unrestricted corrections and reflected mean-correct profiles. It does not cover arbitrary sources or certificate families.',
        'All four pair-transfer moments remain nonzero in every selected current witness. No pair-transfer correction has yet been integrated into the energy certificate.',
        'The averaged signed-charge law has a stationary classical Markov extension. This is not a stationary quantum extension or a globally fixed-particle-number construction.',
        'General quantum representability, molecular or higher-dimensional transfer, and requested-accuracy scalability remain unproved.',
    ]
    summary = {
        'accepted': True, 'cases': cases, 'receipts': receipt_paths,
        'receipt_count': len(receipt_paths), 'source_hash_entries_checked': checks,
        'local_blocks': 94, 'max_local_psd_dimension': 200, 'local_dimension_sum': 4096,
        'total_psd_checks_per_energy': 151, 'sparse_entry_cap': 64, 'family_source_cap': 85,
        'new_independent_spectator_directions': 14, 'combined_hopping_rank': 17,
        'next_pair_transfer_rank': 4,
        'focused_validation': {'tests_passed': 28, 'elapsed_seconds': 111.41},
        'full_validation': validation, 'limitations': limitations,
        'report': 'research/marginal_spectator_hopping.md',
        'numerical_provenance': relative(BASE / 'numerical_provenance.json'),
        'scope': 'Exact open-chain energy intervals, strict escape from the previous fixed spin family, current enlarged-family caps and exact stationary consistency obstructions. Goal remains active.',
    }
    write(BASE / 'numerical_provenance.json', provenance)
    write(BASE / 'combined_summary.json', summary)
    status_path = ROOT / 'results/marginal_final_validation.json'
    status = read(status_path)
    for key in ('command', 'tests_passed', 'subtests_passed', 'elapsed_seconds'):
        status[key] = validation[key]
    status['latest_validation_scope'] = 'Full current suite after ENERGYv13/FAMILYv9 spectator hopping integration; 847 tests and 102 subtests pass, JUnit confirms zero errors/failures/skips. Historical sections retain their original validation evidence.'
    status['full_suite_scope'] = status['latest_validation_scope']
    status['spectator_hopping'] = summary
    write(status_path, status)
    print(json.dumps({'receipts': len(receipt_paths), 'source_hash_entries_checked': checks, 'tests': validation['tests_passed']}))


if __name__ == '__main__':
    main()
