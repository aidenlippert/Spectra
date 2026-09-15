"""Assemble cost/dependency/error records without inventing missing timings."""
from fractions import Fraction as F
import json
from pathlib import Path
from research.global_response_20260913 import global_program as g


def read(name): return json.loads((g.OUT/name).read_text())


def main():
    integrated = read('integrated.json'); discovery = read('global_discovery.json')
    time_rows = []
    for r in discovery['cases']:
        time_rows.append({'stage': r['case']+' uniform response construction',
                          'seconds': r['discovery']['construction_seconds'],
                          'source': 'global_discovery.json', 'includes_shared_input_discovery': False})
    time_rows += [
        {'stage': 'fresh response construction', 'seconds': read('fresh_response_discovery.json')['discovery']['construction_seconds'], 'source': 'fresh_response_discovery.json'},
        {'stage': 'fresh integrals/HF/FCI/tail generation', 'seconds': read('fresh_h6_1p73/construction.json')['wall_seconds'], 'source': 'fresh_h6_1p73/construction.json'},
        {'stage': 'uniform exact replay, first three cases', 'seconds': read('global_replay.json')['wall_seconds'], 'source': 'global_replay.json'},
        {'stage': 'all-block H6 diagnosis and matched execution', 'seconds': read('execution.json')['wall_seconds'], 'source': 'execution.json'},
        {'stage': 'terminal-aware H6 execution plus preparation', 'seconds': read('tuned_execution.json')['wall_seconds'], 'source': 'tuned_execution.json'},
        {'stage': 'independent inherited-SOS terminal/target diagnosis', 'seconds': read('reference_diagnosis.json')['wall_seconds'], 'source': 'reference_diagnosis.json'},
        {'stage': 'final integrated exact construction/replay', 'seconds': integrated['wall_seconds'], 'source': 'integrated.json', 'peak_RSS_bytes': integrated['peak_RSS_bytes']},
    ]
    for name in ('frozen_h6', 'fresh_h6_1p73', 'direct_h6_64', 'direct_h6_256_b2', 'direct_h8_64'):
        path = g.OUT/'joint_factor'/name/'receipt.json'
        if path.exists():
            r = json.loads(path.read_text()); time_rows.append({'stage': name+' factor construction', 'seconds': r['wall_seconds'], 'source': str(path.relative_to(g.OUT))})
    for name in ('h6', 'h8'):
        r = read(f'interference/{name}_weighted.json')
        time_rows.append({'stage': name+' weighted block search', 'seconds': r['seconds'], 'source': f'interference/{name}_weighted.json'})
    for r in read('upper_model/compact_trials.json'):
        time_rows.append({'stage': r['case']+' Slater optimization', 'seconds': r['seconds'], 'objective_evaluations': r['objective_calls'], 'source': 'upper_model/compact_trials.json'})
    for r in read('lifted_upper/discovery.json')['cases']:
        c = r['certificate']; time_rows.append({'stage': 'lift '+c['fixture_sha256'][:12], 'seconds': c['construction_seconds'], 'source': 'lifted_upper/discovery.json'})
    old_backend = read('flatten/benchmark.json')
    for key in ('nested_seconds', 'flattened_seconds'):
        time_rows.append({'stage': 'earlier Python-list/Fraction-affine backend '+key,
                          'seconds': old_backend[key], 'source': 'flatten/benchmark.json',
                          'scope': 'A different implementation from the main vectorized matched benchmark; do not mix speed ratios.'})
    for name, key in (('exceptional/summary.json', 'construction_seconds'), ('exceptional/low_exception_summary.json', 'wall_seconds')):
        for row in read(name)['rows']:
            time_rows.append({'stage': row['case']+' '+name, 'seconds': row[key], 'source': name})
    for name in ('frozen_h6', 'fresh_h6_1p73'):
        sr = read('joint_factor/spectral_local/'+name+'/receipt.json')
        time_rows.append({'stage': name+' post-hoc spectral residual analysis', 'seconds': sr['wall_seconds'],
                          'source': 'joint_factor/spectral_local/'+name+'/receipt.json',
                          'scope': 'Candidate analysis timing; final acceptance expands the unchanged original certificate.'})
    dependency_rows = [
        {'component': 'Hamiltonian fixtures', 'role': 'Ordinary rational one/two-body input data; orbital/integral generation still costs work.', 'many_body_dependency': 'None in data format; FCI reference generation is separate.'},
        {'component': 'Target b', 'role': 'Frozen energy target, selected relative to an enumerated reference upper in tight comparisons.', 'many_body_dependency': 'Inherited scalar selected from FCI upper; it is not a newly discovered compact upper.'},
        {'component': 'Shared density tail', 'role': 'Orbital density factors and exact coefficient residual.', 'many_body_dependency': 'Coefficient-space spectral proposal, not full fixed-N enumeration.'},
        {'component': 'Uniform global response', 'role': 'New gap tangents, local CAR coupling and scalar Chebyshev program.', 'many_body_dependency': 'No full fixed-N enumeration in construction/replay; benchmark execution uses full block vectors.'},
        {'component': 'H6/H8 terminal closure', 'role': 'Inherited full cubic SOS and exterior residual factors, rechecked exactly.', 'many_body_dependency': 'No full fixed-N enumeration, but large original Gram discovery. H6/H8 recorded historical discovery: 259.309/721.563 s.', 'historical_source': 'research/certificate_scaling/CUBIC_PRECISION_RESULTS.md'},
        {'component': 'H6 1.6 Angstrom closure', 'role': 'Inherited expanded retained proof, transported after 0.010 mHa target relaxation.', 'many_body_dependency': '924 fixed-N labels and retained matrices still constructed.'},
        {'component': 'Tight upper witnesses', 'role': '200 amplitudes for H6 fixtures, 2468 for improved H8.', 'many_body_dependency': 'FCI discovery remains. Exact replay streams support times Hamiltonian terms.'},
        {'component': 'Slater upper control', 'role': 'Exact norm/projector/Wick contraction from rational parameters.', 'many_body_dependency': 'No many-body vector on accepting path; small explicit oracle only in tests.'},
        {'component': 'Response lift upper', 'role': 'Single QH excitation and exact moments.', 'many_body_dependency': '140/699 distinct intermediate determinant labels on H6/H8; no full fixed-N basis.'},
        {'component': 'CH2 comparison', 'role': 'Compact pure-spin upper controls versus preserved exact model lower proofs.', 'many_body_dependency': 'The inherited spin lower proofs still enumerate 400/225 spin-projection labels.'},
    ]
    errors = []
    for row in integrated['reference_closed']:
        errors.append({'case': row['case'], 'target_Ha': row['lower_Ha'],
                       'uniform_eta_Ha': row['uniform_response']['residual_penalty_Ha'],
                       'uniform_terminal_margin_Ha': row['uniform_terminal']['terminal_margin_Ha'],
                       'tuned_eta_Ha': row['terminal_aware_response']['residual_penalty_Ha'],
                       'tuned_terminal_margin_Ha': row['terminal_aware_margin']['terminal_margin_Ha'],
                       'error_metric': 'Identity on current retained sector; each nested shift paid once.',
                       'inherited_residual_allowances': 'Already included in the replayed inherited lower; never omitted or counted as zero.'})
    result = {'stage_times': time_rows, 'times_are_not_an_additive_end_to_end_total': True,
              'dependencies': dependency_rows, 'accepted_error_budgets': errors,
              'additional_error_rules': ['Rational orbital PSD and scalar endpoint rounding are outward.',
                'Floating optimization/eigenvalues do not accept lower bounds.',
                'Slater and lifted upper quotients have exact positive norms; no rounding allowance is omitted.',
                'The 1.6 Angstrom transported target pays an explicit 1/100000 Ha loss.',
                'Floating execution is diagnostic; exact acceptance does not rely on its roundoff.',
                'No allowance for basis, geometry, frozen-core or other physical-model error is certified here.'],
              'unavailable_costs_and_rejected_attempts': [
                {'attempt': 'Early scaled-H flatten benchmark', 'status': 'Invalidated; excluded from speedup claims.', 'elapsed_seconds': None, 'evidence': 'flatten/INVALIDATED_SCALED_BENCHMARK.md'},
                {'attempt': 'Initial overly broad weighted search', 'status': 'Stopped; replaced by radius<=1, <=2187 candidates.', 'elapsed_seconds': None},
                {'attempt': 'Repeated-decode upper optimizer', 'status': 'Abandoned after exceeding command window; replaced by predecoded bounded controls.', 'elapsed_seconds': None},
                {'attempt': 'Wrong-map and zero-coupling exceptional claims', 'status': 'Rejected. Corrected local CAR replay and nonzero-transition tests are authoritative.', 'elapsed_seconds': None, 'evidence': 'exceptional/invalidated_zero_coupling/INVALIDATED.md'},
                {'attempt': 'Retagged local residual analysis', 'status': 'Original certificates unchanged. Final replay expands the ORIGINAL certificate and checks body order before accepting candidate spectral factors.', 'elapsed_seconds': None},
              ],
              'campaign_total_runtime': None, 'initial_30_minute_budget_compliance_fully_audited': False,
              'missing_time_is_not_zero': True, 'new_complete_cheap_terminal_discovery': False}
    (g.OUT/'accounting.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'measured_stage_records': len(time_rows), 'total_runtime_claimed': False}))


if __name__ == '__main__': main()
