"""Run declared stages serially, with separate durable attempt costs and limits."""
import argparse
import json
from pathlib import Path
import subprocess
import time
from research.interacting_scaling_20260915.budget import ROOT, OUT, STD, NUM, CHEM, dump

PREFIX = 'research.interacting_scaling_20260915.'


def steps_for(case, discover_state=False, construct_nonsinglet=False, solve_seconds=300, prepared=False, direct_nonsinglet=False):
    if construct_nonsinglet and direct_nonsinglet:
        raise ValueError('Choose one declared magnetic constructor')
    steps = []
    if discover_state:
        command = [NUM, '-B', '-m', PREFIX+'state', str(case)]
        if (case/'continuation_seed.json').exists():
            command += ['--initial', str(case/'continuation_seed.json')]
        steps.append(('state', 900, command))
    if not (case/'upper.json').exists():
        steps.append(('upper', 180, [STD, '-B', '-S', '-m', 'research.transfer_solver_20260915.actions', str(case), 'upper']))
    if construct_nonsinglet:
        steps.append(('nonsinglet', 300, [NUM, '-B', '-m', 'research.transfer_solver_20260915.actions', str(case), 'nonsinglet']))
    if direct_nonsinglet:
        magnetic = OUT/'cases'/(case.name+'_magnetic')
        steps += [('magnetic_initialize', 20, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'initialize', str(case), magnetic.name]),
            ('magnetic_prepare', 300, [NUM, '-B', '-m', PREFIX+'prepare', str(magnetic)]),
            ('magnetic_solve', 195, [NUM, '-B', '-m', PREFIX+'solve', str(magnetic), 'solve', '--seconds', '150']),
            ('magnetic_check', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'check', str(magnetic), 'solve', str(magnetic/'exact.json')]),
            ('magnetic_attach', 180, [STD, '-B', '-S', '-m', PREFIX+'nonsinglet', 'attach', str(case), str(magnetic)])]
    if not prepared:
        steps.append(('prepare', 900, [NUM, '-B', '-m', PREFIX+'prepare', str(case)]))
    steps.append(('solve', min(900, solve_seconds+45), [NUM, '-B', '-m', PREFIX+'solve', str(case), 'solve', '--seconds', str(solve_seconds)]))
    steps.append(('replay', 600, [STD, '-B', '-S', '-m', 'research.nvidia_followup_20260915.strict_replay', str(case), 'solve', str(case/'exact')]))
    return steps


def execute(name, steps, deadline=None):
    destination = OUT/'executions'/name
    destination.mkdir(parents=True, exist_ok=False)
    dump(destination/'protocol.json', {'steps': steps, 'one_heavy_process_at_a_time': True,
        'host': 'local Mac; one BLAS/Numba thread', 'external_spending': False})
    start = time.monotonic()
    outcomes = []
    for label, seconds, command in steps:
        if deadline is not None:
            remaining = int(deadline-time.monotonic())
            if remaining < 1:
                outcomes.append({'stage': label, 'exit_code': 124, 'not_started': 'Complete campaign budget exhausted'})
                break
            seconds = min(seconds, remaining)
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', name+'_'+label,
            '--seconds', str(seconds), '--']+command, cwd=ROOT)
        outcomes.append({'stage': label, 'exit_code': result.returncode, 'effective_cap_seconds': seconds})
        if result.returncode:
            break
    record = {'steps': outcomes, 'wall_seconds': time.monotonic()-start,
        'completed': len(outcomes) == len(steps) and not outcomes[-1]['exit_code'],
        'last_stage': outcomes[-1]['stage']}
    dump(destination/'execution.json', record)
    print(json.dumps(record), flush=True)
    return record


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('name')
    p.add_argument('case', type=Path)
    p.add_argument('--state', action='store_true')
    p.add_argument('--nonsinglet', action='store_true')
    p.add_argument('--direct-nonsinglet', action='store_true')
    p.add_argument('--solve-seconds', type=int, default=300)
    a = p.parse_args()
    execute(a.name, steps_for(a.case.resolve(), a.state, a.nonsinglet, a.solve_seconds, direct_nonsinglet=a.direct_nonsinglet))
