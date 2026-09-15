"""Matched magnetic-component control for the same spin-averaged cone.

Only highest-weight reduction is disabled. The input-derived dictionaries,
coefficient equations, cold primal point, MPS functional and solve rule remain.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, PREFIX


def prepare_control(case):
    import numpy as np
    from research.transfer_solver_20260915 import prepare
    prepare.highest_weights = lambda groups: [(i, np.eye(len(g['words'])),
                                               {'kind': 'unreduced_magnetic_component'}) for i, g in enumerate(groups)]
    prepare.build(case)


def run(source):
    started = time.monotonic()
    case = OUT/'controls'/(source.name+'_magnetic')
    case.mkdir(parents=True, exist_ok=False)
    (case/'mps').mkdir()
    inputs = ['fixture.json', 'upper.json', 'nonsinglet.json', 'mps/state.json']
    hashes = {}
    for name in inputs:
        shutil.copyfile(source/name, case/name)
        hashes[name] = hashlib.sha256((case/name).read_bytes()).hexdigest()
    common = lambda module, args: [NUM, '-B', '-m', PREFIX+module]+args
    steps = [('prepare', 900, common('representation_control', ['prepare', str(case)])),
             ('stage1', 400, common('solve', [str(case), 'stage1', '--seconds', '300', '--mu', '2'])),
             ('stage2', 300, common('solve', [str(case), 'stage2', '--seconds', '180', '--mu', '.03',
                                           '--restart', str(case/'stage1/checkpoint.npz')])),
             ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'actions', str(case), 'replay'])]
    dump(case/'protocol.json', {'source_case': str(source), 'input_hashes': hashes,
                                'identical_cold_physical_initial_primal': 'zero Gram and ideal coordinates',
                                'attainable_singlet_cone': 'Same spin-averaged mixed dictionary; highest-weight decomposition disabled',
                                'common_state_and_nonsinglet_costs_excluded_from_both_representation_comparisons': True,
                                'steps': [{'name': n, 'seconds': t, 'command': c} for n, t, c in steps]})
    status = 'completed'
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', case.name+'_'+name,
                                 '--seconds', str(seconds), '--']+command, cwd=ROOT)
        if result.returncode:
            status = 'stopped_after_failed_step'
            break
    dump(case/'execution.json', {'status': status, 'seconds': time.monotonic()-started, 'last_step': name,
                                'shared_dependency_costs_excluded': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'prepare'))
    p.add_argument('case', type=Path)
    a = p.parse_args()
    prepare_control(a.case.resolve()) if a.action == 'prepare' else run(a.case.resolve())
