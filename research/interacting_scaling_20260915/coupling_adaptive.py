"""Apply the frozen nested proof rule along a globally charge-aware coupling path."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.pipeline import execute, steps_for, PREFIX
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, dump
from research.interacting_scaling_20260915.nonsinglet import construct as screen


def run(source, name, original, total_seconds=6000):
    frozen = json.loads((OUT/'heldout_protocol.json').read_text())
    directory = OUT/'coupling'/name
    directory.mkdir(parents=True, exist_ok=False)
    def require_frozen():
        for file, value in frozen['source_sha256'].items():
            if hashlib.sha256((ROOT/'research/interacting_scaling_20260915'/file).read_bytes()).hexdigest() != value:
                raise ValueError('The frozen adaptation procedure changed')
    require_frozen()
    dump(directory/'protocol.json', {'source': str(source), 'original': str(original),
        'couplings': ['0', '1/4', '1/2', '1'], 'fragment_width': 2,
        'levels': frozen['levels'], 'per_level_solve_seconds': frozen['per_level_solve_seconds'],
        'total_continuation_budget_seconds': total_seconds,
        'global_particle_number_only': True, 'coupling_changes_charge_transfer_terms': True,
        'source_integral_and_rotation_costs_additional': True,
        'fresh_state_at_zero_then_continuation': True, 'earlier_Gram_solutions_not_transferred_between_couplings': True,
        'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'frozen_procedure_sha256': hashlib.sha256((OUT/'heldout_protocol.json').read_bytes()).hexdigest()})
    start = time.monotonic(); deadline = start+total_seconds
    previous_state, outcomes = None, []
    for coupling, label in [('0', 'l0'), ('1/4', 'l025'), ('1/2', 'l05'), ('1', 'l1')]:
        base = OUT/'cases'/(name+'_'+label+'_base')
        command = [STD, '-B', '-S', '-m', PREFIX+'coupling', str(source), base.name, coupling]
        if previous_state: command += ['--previous', str(previous_state)]
        initial = execute(base.name+'_init', [('initialize', 20, command)], deadline)
        row = {'lambda': coupling, 'base_case': str(base), 'initialization': initial, 'levels': []}
        outcomes.append(row)
        if not initial['completed']: break
        common = steps_for(base, True, direct_nonsinglet=True)
        common = common[:next(i for i, (label, _, _) in enumerate(common) if label == 'magnetic_initialize')]
        state = execute(base.name+'_state', common, deadline)
        row['state_and_upper'] = state
        if not state['completed']: break
        previous_state = base
        magnetic = screen(base, deadline); row['screen'] = magnetic
        if not magnetic['completed']: continue
        previous_proof, best = None, None
        for level, widths in enumerate(frozen['levels']):
            require_frozen()
            if max(widths) > json.loads((base/'fixture.json').read_text())['modes']//2: break
            case = OUT/'cases'/f'{name}_{label}_level{level}'
            command = [STD, '-B', '-S', '-m', PREFIX+'cases', case.name, str(base), '--widths']+list(map(str, widths))
            command += ['--collective-widths']+list(map(str, [k for k in widths if k >= 3]))
            stages = [('initialize', 20, command), ('prepare', 900, [NUM, '-B', '-m', PREFIX+'prepare', str(case)])]
            if previous_proof:
                stages.append(('transport', 180, [NUM, '-B', '-m', PREFIX+'warmstart', str(previous_proof), 'solve', str(case)]))
            seconds = frozen['per_level_solve_seconds']
            command = [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', str(seconds)]
            if previous_proof: command += ['--restart', str(case/'nested_start.npz')]
            stages.append(('solve', min(900, seconds+45), command))
            if coupling == '1':
                command = [STD, '-B', '-S', '-m', PREFIX+'complete', str(original), str(case), 'solve', str(case/'exact')]
                interval_path = case/'original_interval.json'
            else:
                command = [STD, '-B', '-S', '-m', 'research.nvidia_followup_20260915.strict_replay', str(case), 'solve', str(case/'exact')]
                interval_path = case/'exact/interval.json'
            stages.append(('replay', 600, command))
            result = execute(case.name, stages, deadline)
            attempt = {'level': level, 'case': str(case), 'execution': result}
            row['levels'].append(attempt)
            if not result['completed']: break
            interval = json.loads(interval_path.read_text())
            attempt['interval'] = interval
            if best is None or interval['width_mHa'] < best['interval']['width_mHa']: best = attempt
            previous_proof = case
            if best['interval']['target_met']: break
        row['best'] = best
        row['target_met'] = bool(best and best['interval']['target_met'])
    dump(directory/'results.json', {'outcomes': outcomes, 'wall_seconds': time.monotonic()-start,
        'fully_coupled_target_met': any(d['lambda'] == '1' and d.get('target_met') for d in outcomes),
        'all_state_and_failed_proof_continuation_costs_included': True,
        'original_integrals_and_orbitals_are_additional': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('source', type=Path); p.add_argument('name'); p.add_argument('original', type=Path)
    a = p.parse_args(); run(a.source.resolve(), a.name, a.original.resolve())
