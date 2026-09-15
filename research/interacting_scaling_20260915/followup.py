"""Bounded remaining validation jobs; run only after the matched cold comparisons."""
import json
from pathlib import Path
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.pipeline import execute, PREFIX


def run():
    for label in ('h8_matched_direct', 'h8_matched_reference', 'h10_matched_direct', 'h10_matched_reference'):
        if not (OUT/'cold'/label/'result.json').exists():
            raise ValueError('Finish the matched cold comparisons before this queue')
    first = execute('fragment_control_tests', [('unit', 60, [STD, '-B', '-S', '-m', 'unittest', PREFIX+'test_fragment_control', PREFIX+'test_capacity'])])
    if not first['completed']: raise RuntimeError('Resolve the new control/envelope test failures before validation')
    execute('immutable_bundle_replay', [('h4', 60, [STD, '-B', '-S', '-m', PREFIX+'replay_bundle',
        str(OUT/'models/h4_cold_integration'), str(OUT/'cases/h4_cold_integration_level0'), str(OUT/'independent_bundle_replay/h4')])])
    case, folder = OUT/'cases/h4_coupling_l0', OUT/'fragment_controls/h4'
    execute('h4_charge_complete_control', [('construct', 30, [NUM, '-B', '-m', PREFIX+'fragment_control', 'construct', str(case), str(folder)]),
        ('check', 30, [STD, '-B', '-S', '-m', PREFIX+'fragment_control', 'check', str(case), str(folder)])])
    from research.interacting_scaling_20260915.heldout import run as heldout
    heldout('h12_heldout')
    from research.interacting_scaling_20260915.capacity import run_conditionally
    h12 = run_conditionally()
    for name in ('water_heldout', 'h4_631g_heldout'): heldout(name)
    if h12['target_met']:
        from research.interacting_scaling_20260915.cold import run as cold
        reference_budget = json.loads((OUT/'reference_budget_protocol.json').read_text())
        frozen = json.loads((OUT/'heldout_protocol.json').read_text())
        cold('h12_heldout_reference', frozen['models']['h12_heldout'],
            solve_seconds=reference_budget['main_solver_cap_seconds'],
            total_seconds=reference_budget['complete_pipeline_cap_seconds'], reference_mode=True)
    # The first gate is measured, not inferred from the combinatorial state count.
    gate = {'H12_target_met': h12['target_met'], 'H12_complete_seconds': h12['wall_seconds'],
        'H16_attempted': False, 'required': 'H12 target within 1800 s, then a representation/memory forecast within the frozen local envelope',
        'reason': 'H12 did not meet the target within 1800 seconds' if not h12['target_met'] or h12['wall_seconds'] > 1800 else 'Structural forecast required before any H16 allocation'}
    dump(OUT/'h16_gate.json', gate)
    from research.interacting_scaling_20260915.coupling_adaptive import run as coupling
    coupling(OUT/'cases/h8_matched_direct_base', 'h8_adaptive_coupling', OUT/'models/h8_matched_direct')
    case, folder = OUT/'cases/h8_adaptive_coupling_l0_base', OUT/'fragment_controls/h8'
    execute('h8_charge_complete_control', [('construct', 30, [NUM, '-B', '-m', PREFIX+'fragment_control', 'construct', str(case), str(folder)]),
        ('check', 30, [STD, '-B', '-S', '-m', PREFIX+'fragment_control', 'check', str(case), str(folder)])])
    from research.interacting_scaling_20260915.baselines import numerical, exact
    for name in ('h8_matched_direct', 'h10_matched_direct', 'h12_heldout', 'water_heldout', 'h4_631g_heldout'):
        numerical(OUT/'models'/name)
    for name in ('h4_cold_integration', 'water_heldout', 'h4_631g_heldout'):
        exact(OUT/'models'/name, name+'_enumerated_control')
    tests = ['test_coordinates', 'test_dictionary', 'test_singlet_trace', 'test_spin_patterns',
        'test_dual_repair', 'test_nonsinglet', 'test_pipeline', 'test_fragment_control', 'test_capacity']
    execute('final_regressions', [('unit', 120, [NUM, '-B', '-m', 'unittest']+[PREFIX+t for t in tests])])
    print('All queued validation attempts finished; review exact outcomes and the H16 gate before accounting and sealing.', flush=True)


if __name__ == '__main__': run()
