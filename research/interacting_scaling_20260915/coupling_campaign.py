"""Bounded continuation through exact interfragment coupling, including lambda=1."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from research.interacting_scaling_20260915.pipeline import execute, steps_for, PREFIX
from research.interacting_scaling_20260915.budget import OUT, STD, dump


def run(source, name, original=None):
    directory = OUT/'coupling'/name
    directory.mkdir(parents=True, exist_ok=False)
    dump(directory/'protocol.json', {'source': str(source), 'couplings': ['0', '1/4', '1/2', '1'],
        'fragment_width': 2, 'local_widths': [2, 3], 'collective_pair_supports': True,
        'state': 'Fresh global-N Neel start at zero; explicitly rebound continuation states thereafter',
        'nonsinglet': 'Direct quadratic, exact occupation-trace initializer',
        'per_singlet_solve_seconds': 120, 'global_full_cubic_map': False})
    start = time.monotonic()
    previous = None
    outcomes = []
    for coupling, label in [('0', 'l0'), ('1/4', 'l025'), ('1/2', 'l05'), ('1', 'l1')]:
        case = OUT/'cases'/(name+'_'+label)
        command = [STD, '-B', '-S', '-m', PREFIX+'coupling', str(source), case.name, coupling]
        if previous: command += ['--previous', str(previous)]
        init = execute(case.name+'_init', [('initialize', 20, command)])
        if not init['completed']: raise RuntimeError('Coupling construction did not finish')
        result = execute(case.name, steps_for(case, True, False, 120, direct_nonsinglet=True))
        record = {'lambda': coupling, 'case': str(case), 'pipeline': result}
        if result['completed'] and (case/'exact/complete_replay.json').exists():
            if not json.loads((case/'exact/complete_replay.json').read_text())['all_upper_and_lower_dependencies_rechecked']:
                raise ValueError('Complete coupling replay is required')
            interval = json.loads((case/'exact/interval.json').read_text())
            record.update(width_mHa=interval['width_mHa'], target_met=interval['target_met'])
        if result['completed'] and coupling == '1' and original:
            transfer = execute(case.name+'_original', [('transfer', 600, [STD, '-B', '-S', '-m', PREFIX+'rotation',
                'transfer', str(original), str(case), str(case/'exact'), str(case/'original_interval.json')])])
            record['original_transfer'] = transfer
            if transfer['completed']:
                record['original_interval'] = json.loads((case/'original_interval.json').read_text())
        outcomes.append(record)
        if (case/'mps/state.json').exists(): previous = case
    dump(directory/'results.json', {'outcomes': outcomes, 'wall_seconds': time.monotonic()-start,
        'continuation_costs_all_included': True, 'integrals_and_source_rotation_are_additional': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('source', type=Path); p.add_argument('name'); p.add_argument('--original', type=Path)
    a = p.parse_args(); run(a.source.resolve(), a.name, a.original.resolve() if a.original else None)
