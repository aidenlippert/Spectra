"""Cold HF-product initialization control for the observed random-state bottleneck."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, PREFIX, CASES
from research.molecular_collective_20260913.core import digest


def state(case):
    data = json.loads((case/'fixture.json').read_text())
    counts = [data['particles']//2]*2
    charges, tensors = [[[0, 0]]], []
    current = [0, 0]
    for i in range(data['modes']):
        bit = int(i//2 < counts[i % 2])
        tensors.append([[0, bit, 0, 1]])
        current = current.copy()
        current[i % 2] += bit
        charges.append([current])
    cert = {'kind': 'integer_charge_mps_v1', 'fixture_sha256': digest(data),
        'modes': data['modes'], 'particles': data['particles'], 'spin_counts': counts,
        'denominator': 1, 'bond_charges': charges, 'tensors': tensors}
    from research.correlated_pair_20260913.mps_exact import State
    State(data, cert)
    seed = case/'hf_seed.json'
    dump(seed, cert)
    from research.correlated_pair_20260913.mps_spatial import run
    folder = case/'mps'
    folder.mkdir(exist_ok=False)
    run(data, counts, folder, bond=64, sweeps=6, initial=str(seed))
    dump(folder/'initialization_clarification.json', {'actual_initialization': 'New exact HF product from declared RHF orbital occupations',
        'old_many_body_state_used': False, 'initial_path_is_this_runs_fresh_one_determinant_state': True,
        'note': 'The imported proposer has a generic random-initializer description; its initial_MPS_path and this record specify the actual supplied initializer.'})


def run():
    source = OUT/'adaptive/h10_bond64'
    old = json.loads((source/'protocol.json').read_text())
    old_execution = json.loads((source/'execution.json').read_text())
    if old_execution['status'] == 'completed':
        raise ValueError('This extra initialization control is reserved for an unfinished size retry')
    case = OUT/'adaptive/h10_hf64'
    case.mkdir(parents=True, exist_ok=False)
    dump(case/'specification.json', {**CASES['h10_size'], 'role': 'HF product initialization after measured random-state construction bottleneck'})
    steps = []
    for row in old['steps']:
        command = [x.replace(str(source), str(case)) for x in row['command']]
        if row['name'] == 'state':
            command = [NUM, '-B', '-m', 'research.transfer_followup_20260915.hf_size', 'state', str(case)]
        steps.append({**row, 'command': command})
    dump(case/'protocol.json', {'preceding_attempt': str(source), 'preceding_execution': old_execution,
        'new_integrals_new_HF_product_no_old_state': True, 'bond_cap': 64, 'sweeps_cap': 6,
        'same_lower_family_optimizer_and_accepting_checker': True, 'steps': steps,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'not_a_success_of_the_original_random_state_frozen_rule': True})
    env = os.environ.copy()
    env['NUMBA_CACHE_DIR'] = str(case/'numba_cache')
    start = time.monotonic()
    status = 'completed'
    for row in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h10_hf64_'+row['name'],
                                 '--seconds', str(row['seconds']), '--']+row['command'], cwd=ROOT, env=env)
        if result.returncode:
            status = 'stopped_after_failed_step'
            break
    dump(case/'execution.json', {'status': status, 'last_step': row['name'], 'elapsed_seconds': time.monotonic()-start,
                                 'earlier_failed_attempt_costs_are_additional': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'state'))
    p.add_argument('case', type=Path, nargs='?')
    a = p.parse_args()
    run() if a.action == 'run' else state(a.case.resolve())
