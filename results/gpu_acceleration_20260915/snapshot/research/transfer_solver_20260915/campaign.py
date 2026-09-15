"""Account for a complete cold construction on a predeclared molecular input."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump

STD = '/opt/homebrew/Caskroom/miniconda/base/bin/python'
NUM = str(ROOT/'.venv-correlated/bin/python')
CHEM = str(ROOT/'.venv-molecule/bin/python')
PREFIX = 'research.transfer_solver_20260915.'
CASES = {
    'h4_control': {'geometry': [['H', [0., 0., 1.4*i]] for i in range(4)], 'basis': 'sto-3g', 'role': 'development control'},
    'h8_cold': {'geometry': [['H', [0., 0., 1.4*i]] for i in range(8)], 'basis': 'sto-3g', 'role': 'cold reconstruction of known geometry'},
    'h6_asymmetric': {'geometry': [['H', [0.07*i*i, 0.03*(i % 3), 1.4*i]] for i in range(6)], 'basis': 'sto-3g', 'role': 'symmetry-breaking transfer'},
    'water_asymmetric': {'geometry': [['O', [0., 0., 0.]], ['H', [0.7586, 0., 0.5043]], ['H', [-0.80, 0.02, 0.53]]],
                         'basis': 'sto-3g', 'role': 'different molecule and electron/orbital counts'},
    'h10_size': {'geometry': [['H', [0., 0., 1.4*i]] for i in range(10)], 'basis': 'sto-3g', 'role': 'modest size increase'},
}


def plan(case):
    module = lambda py, name, args: [py, '-B']+(['-S'] if py == STD else [])+['-m', PREFIX+name]+args
    c = str(case)
    return [('generate', 90, module(CHEM, 'generate', [c])),
            ('state', 900, module(NUM, 'actions', [c, 'state'])),
            ('upper', 600, module(STD, 'actions', [c, 'upper'])),
            ('nonsinglet', 240, module(NUM, 'actions', [c, 'nonsinglet'])),
            ('prepare', 900, module(NUM, 'prepare', [c])),
            ('stage1', 400, module(NUM, 'solve', [c, 'stage1', '--seconds', '300', '--mu', '2'])),
            ('stage2', 300, module(NUM, 'solve', [c, 'stage2', '--seconds', '180', '--mu', '.03',
                                             '--restart', str(case/'stage1/checkpoint.npz')])),
            ('replay', 600, module(STD, 'actions', [c, 'replay']))]


def run(name):
    start = time.monotonic()
    case = OUT/'cases'/name
    case.mkdir(parents=True, exist_ok=False)
    dump(case/'specification.json', CASES[name])
    sources = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted((ROOT/'research/transfer_solver_20260915').glob('*.py'))}
    steps = [{'name': n, 'timeout_seconds': t, 'command': c} for n, t, c in plan(case)]
    dump(case/'protocol.json', {'created_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                               'steps': steps, 'sources': sources, 'target_total_width_Ha': '1/625',
                               'old_results_used_as_constructor_inputs': False, 'old_states_or_checkpoints_used': False,
                               'model': CASES[name], 'initial_primal': 'zero',
                               'singlet_method': 'full input-derived mixed highest weights, all mixed linear words, full other dictionaries',
                               'nonsinglet_rank': 4, 'MPS_sweeps': 10, 'MPS_bond_rule': 'min(144,4**floor(spatial_orbitals/2))'})
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['NUMBA_CACHE_DIR'] = str(case/'numba_cache')
    status = 'completed'
    for row in steps:
        for p, digest in sources.items():
            if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != digest:
                raise ValueError('Source changed after this case was frozen')
        cmd = [STD, '-B', '-S', '-m', PREFIX+'budget', '--name', name+'_'+row['name'],
               '--seconds', str(row['timeout_seconds']), '--']+row['command']
        result = subprocess.run(cmd, cwd=ROOT, env=env)
        if result.returncode:
            status = 'stopped_after_failed_step'
            break
    dump(case/'execution.json', {'status': status, 'last_step': row['name'],
                                'elapsed_seconds': time.monotonic()-start,
                                'from_new_integrals': True, 'historical_artifacts_loaded': False,
                                'development_and_failed_other_case_costs_are_separate': True})
    return 0 if status == 'completed' else 1


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', choices=CASES)
    raise SystemExit(run(p.parse_args().case))
