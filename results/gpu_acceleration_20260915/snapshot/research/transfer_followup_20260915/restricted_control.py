"""Declared nested-family control; finite-search failure is not an obstruction."""
import argparse
import hashlib
import shutil
import subprocess
import time
from pathlib import Path
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, PREFIX


def prepare(case):
    from research.transfer_solver_20260915 import prepare as module
    original = module.highest_weights

    def restricted(groups):
        result = []
        for i, basis, receipt in original(groups):
            if receipt['kind'] != 'unchanged_full_dictionary':
                receipt = {**receipt, 'restriction': 'First at most 64 columns of the deterministic rational highest-weight basis',
                           'full_columns': basis.shape[1]}
                basis = basis[:, :64]
            result.append((i, basis, receipt))
        return result

    module.highest_weights = restricted
    module.build(case)


def run():
    started = time.monotonic()
    source = OUT/'cases/h8_cold'
    case = OUT/'controls/h8_restricted64'
    case.mkdir(parents=True, exist_ok=False)
    (case/'mps').mkdir()
    hashes = {}
    for name in ('fixture.json', 'upper.json', 'nonsinglet.json', 'mps/state.json'):
        shutil.copyfile(source/name, case/name)
        hashes[name] = hashlib.sha256((case/name).read_bytes()).hexdigest()
    steps = [('prepare', 900, [NUM, '-B', '-m', 'research.transfer_followup_20260915.restricted_control', 'prepare', str(case)]),
             ('stage1', 400, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'stage1', '--seconds', '300', '--mu', '2']),
             ('stage2', 300, [NUM, '-B', '-m', PREFIX+'solve', str(case), 'stage2', '--seconds', '180', '--mu', '.03',
                               '--restart', str(case/'stage1/checkpoint.npz')]),
             ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'actions', str(case), 'replay'])]
    dump(case/'protocol.json', {'input_hashes': hashes, 'nested_restriction': 'At most first 64 rational columns per mixed highest-weight block',
         'common_initial_primal': 'zero', 'same_state_and_nonsinglet_as_full_control': True,
         'same_optimizer_and_stopping_rule': True, 'different_attainable_family': True,
         'not_claimed_to_be_best_restricted_basis': True, 'no_verified_family_ceiling': True,
         'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'steps': [{'name': n, 'seconds': t, 'command': c} for n, t, c in steps]})
    status = 'completed'
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h8_restricted64_'+name,
                                  '--seconds', str(seconds), '--']+command, cwd=ROOT)
        if result.returncode:
            status = 'stopped_after_failed_step'
            break
    dump(case/'execution.json', {'status': status, 'last_step': name, 'seconds': time.monotonic()-started,
                                 'shared_input_construction_costs_excluded': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'prepare'))
    p.add_argument('case', type=Path, nargs='?')
    a = p.parse_args()
    run() if a.action == 'run' else prepare(a.case.resolve())
