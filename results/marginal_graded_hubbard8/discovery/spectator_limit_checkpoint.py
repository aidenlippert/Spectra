"""Audit the completed limit-refinement receipts and collect their provenance."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'results/marginal_graded_hubbard8/spectator_hopping'
CASES = ('W_zero', 'W_plus_1')
NAMES = ('range_two_replay', 'range_two_family_limit_replay', 'family_escape',
         'spectator_ablation', 'pair_transfer_overlap', 'one_spectator_overlap',
         'spin_overlap', 'coherent_overlap', 'charge_markov_extension',
         'solver_limit_refutation')


def read(path):
    return json.loads(path.read_text())


def relative(path):
    return str(path.relative_to(ROOT))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def main():
    old = read(BASE / 'combined_summary.json')
    receipts = [ROOT / p for p in old['receipts']]
    new_receipts = [BASE / c / 'limit_refined' / (name + '.json') for c in CASES for name in NAMES]
    supplements = [BASE / c / stage / 'range_two_replay.json' for c in CASES for stage in ('cutting', 'subspace48')]
    supplements.append(BASE / 'W_zero/subspace48/range_two_family_limit_replay.json')
    all_receipts = receipts + new_receipts + supplements
    count = 0
    for path in all_receipts:
        data = read(path)
        assert data['accepted'], path
        for name, expected in data['source_sha256'].items():
            assert sha(ROOT / name) == expected, (path, name)
            count += 1
    cases = {}
    for case in CASES:
        folder = BASE / case / 'limit_refined'
        energy = read(folder / 'range_two_replay.json')
        family = read(folder / 'range_two_family_limit_replay.json')
        pair = read(folder / 'pair_transfer_overlap.json')
        gain = F(energy['lower_per_site']) - F(old['cases'][case]['lower_per_site'])
        assert gain > 0
        assert all(p['violated'] for p in pair['separators'])
        assert all(F(v) == 0 for v in read(folder / 'one_spectator_overlap.json')['spectator_moments'].values())
        assert F(energy['upper_per_site']) == F(old['cases'][case]['upper_per_site'])
        cases[case] = {
            'directory': relative(folder), 'target': energy['target'],
            'lower_per_site': energy['lower_per_site'], 'lower_float': float(F(energy['lower_per_site'])),
            'upper_per_site': energy['upper_per_site'], 'upper_float': float(F(energy['upper_per_site'])),
            'lower_gain': str(gain), 'lower_gain_float': float(gain),
            'periodic_lower': family['accepted_periodic_lower'],
            'periodic_family_ceiling': family['periodic_family_upper'],
            'family_gap': family['family_gap'], 'family_gap_float': family['family_gap_float'],
            'previous_family_gap_float': old['cases'][case]['family_gap_float'],
            'family_sources': len(read(folder / 'range_two_family_limit_certificate.json')['mixture']),
            'pair_transfer_separators': pair['separators'],
            'strict_escape_over_spin_family': read(folder / 'family_escape.json')['strict_escape'],
            'invalid_solver_values': read(folder / 'solver_limit_refutation.json')['contradictions'],
            'lower_origin': relative(BASE / case / 'subspace48'),
            'family_origin': relative(BASE / case / ('subspace48' if case == 'W_zero' else 'polished')),
        }
    stages = ('cutting', 'subspace', 'subspace48', 'subspace_atoms', 'atom_pricing',
              'global_subspace', 'global_diagonal', 'scaled_basis', 'slack_repair')
    files = set()
    for case in CASES:
        for stage in stages:
            folder = BASE / case / stage
            if folder.exists():
                files.update(p for p in folder.iterdir() if p.is_file())
        files.update(p for p in (BASE / case).glob('*numeric*.log'))
    drivers = ('spectator_cutting_numeric.py', 'spectator_subspace_numeric.py',
               'spectator_subspace48_numeric.py', 'spectator_subspace_atoms.py',
               'spectator_atom_family_numeric.py', 'spectator_atom_moments.py',
               'spectator_global_subspace_numeric.py', 'spectator_global_diagonal_numeric.py',
               'spectator_scaled_basis_numeric.py', 'spectator_slack_basis_repair.py',
               'spectator_family_numeric.py', 'spectator_solver_limit_refutation.py')
    files.update(ROOT / 'results/marginal_graded_hubbard8/discovery' / name for name in drivers)
    files.add(BASE / 'pre_spectator_sources/subspace48_before_serialization_fix.py')
    provenance = {
        'accepted_mathematical_claim': False,
        'sha256': {relative(p): sha(p) for p in sorted(files)},
        'selected_lower': 'subspace48: 12 rounds, maximum cone dimension48, at most32 active blocks, temporary coefficient radius .05, 15 solver seconds/60 iterations per round. Both reached the round limit, then passed independent exact replay.',
        'selected_W_zero_family': 'spectator_atom_family_numeric.py: subspace48 seed, eight dominant frozen-conic dual vectors and 200 total added rounded/perturbed vectors; scale .0002, determinant anchors, zero pricing rounds. Exact77-source mixture accepted.',
        'selected_W_plus_1_family': 'Unchanged prior polished85-source proposal, copied and freshly matched to the new lower. Historical numerical origin remains in numerical_provenance.json.',
        'failed_pricing': '40 rounds/two eigenvectors with the conic atoms produced83/84-column inconsistent exact bases. Row rescaling100,10000,1000000 did not repair either. Explicit omitted-slack trials produced dependence, negative weights or inconsistency; none accepted.',
        'reuse': 'scaled_basis/candidate_ledger.json.gz retains each bounded finite candidate set for further basis work without rerunning spectral pricing.',
        'scope': 'Provenance only. Numerical status and objective values do not accept physical positivity or family optimality.',
    }
    write(BASE / 'limit_refinement_provenance.json', provenance)
    summary = {
        'accepted': True, 'cases': cases,
        'new_receipts': list(map(relative, new_receipts)),
        'supplemental_receipts': list(map(relative, supplements)),
        'total_receipts_audited_including_previous': len(all_receipts),
        'source_hash_entries_checked': count,
        'production_changed': False,
        'prior_full_regression': old['full_validation'],
        'verification_scope': 'Fresh exact lower/family/obstruction replays and hash audit. No production verifier changed; the prior847-test/102-subtest full regression was not rerun in this refinement.',
        'report': 'research/marginal_spectator_limit_refinement.md',
        'provenance': relative(BASE / 'limit_refinement_provenance.json'),
        'remaining': [
            'The unrestricted current family optimum remains unresolved, especially W=1.',
            'Temporary-box conic duals have nonzero full-family moment residuals. Incomplete global conic objective values are explicitly contradicted by accepted lower certificates.',
            'All four pair-transfer obstructions persist in both selected witnesses; pair-transfer energy integration remains undone.',
            'W=+/-0.1 retain their previous exact checkpoints; no new molecular, long-range or higher-dimensional transfer was established.',
            'General representability and requested-accuracy scalability remain unproved. Goal active.',
        ],
    }
    write(BASE / 'limit_refinement_summary.json', summary)
    status_path = ROOT / 'results/marginal_final_validation.json'
    status = read(status_path)
    status['spectator_limit_refinement'] = summary
    status['latest_validation_scope'] = summary['verification_scope']
    write(status_path, status)
    print(json.dumps({'receipts_audited': len(all_receipts), 'source_hash_entries_checked': count,
                      'cases': {k: {'lower': v['lower_float'], 'family_gap': v['family_gap_float'], 'sources': v['family_sources']} for k,v in cases.items()}}))


if __name__ == '__main__':
    main()
