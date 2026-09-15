"""Canonical lower certificate paired with the independently verified rotated upper.

A newly built Hartree--Fock product supplies proposal moments only. The accepted
upper belongs to the original Hamiltonian through the exact orbital-rotation
receipt; no full-space reference enters either constructor.
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
from research.molecular_collective_20260913.core import digest

SOURCE = OUT/'cases/h10_size'
UPPER = OUT/'rotated_h10'
CASE = OUT/'adaptive/h10_rotated_upper'


def initialize():
    from research.correlated_pair_20260913.mps_exact import State
    data = json.loads((SOURCE/'fixture.json').read_text())
    upper = json.loads((UPPER/'original_upper.json').read_text())
    if upper['original_fixture_sha256'] != digest(data) or not upper['exact_orbital_orthogonality']:
        raise ValueError('Exact original-model upper binding required')
    shutil.copyfile(SOURCE/'fixture.json', CASE/'fixture.json')
    shutil.copyfile(UPPER/'original_upper.json', CASE/'upper.json')
    counts = [data['particles']//2]*2
    current, charges, tensors = [0, 0], [[[0, 0]]], []
    for i in range(data['modes']):
        bit = int(i//2 < counts[i % 2])
        tensors.append([[0, bit, 0, 1]])
        current = current.copy()
        current[i % 2] += bit
        charges.append([current])
    state = {'kind': 'integer_charge_mps_v1', 'fixture_sha256': digest(data),
        'modes': data['modes'], 'particles': data['particles'], 'spin_counts': counts,
        'denominator': 1, 'bond_charges': charges, 'tensors': tensors}
    State(data, state)
    dump(CASE/'mps/state.json', state)
    dump(CASE/'input_roles.json', {
        'proposal_MPS': 'Fresh canonical HF product; used only for moment initialization',
        'upper_proof': str(UPPER/'original_upper.json'),
        'upper_proof_file_sha256': hashlib.sha256((UPPER/'original_upper.json').read_bytes()).hexdigest(),
        'upper_is_NOT_the_HF_guide_energy': True,
        'integrals_and_rotation_state_costs_are_additional': True,
        'FCI_used': False, 'many_body_reference_checkpoint_used': False})


def run():
    if not (UPPER/'original_upper.json').exists():
        raise ValueError('Finish the exact orbital-rotation upper first')
    CASE.mkdir(parents=True, exist_ok=False)
    steps = [
        ('initialize', 30, [STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.rotated_lower', 'initialize']),
        ('nonsinglet', 240, [NUM, '-B', '-m', PREFIX+'actions', str(CASE), 'nonsinglet']),
        ('prepare', 900, [NUM, '-B', '-m', 'research.transfer_followup_20260915.adaptive_size', 'prepare', str(CASE)]),
        ('stage1', 400, [NUM, '-B', '-m', PREFIX+'solve', str(CASE), 'stage1', '--seconds', '300', '--mu', '2']),
        ('stage2', 300, [NUM, '-B', '-m', PREFIX+'solve', str(CASE), 'stage2', '--seconds', '180', '--mu', '.03',
                        '--restart', str(CASE/'stage1/checkpoint.npz')]),
        ('replay', 600, [STD, '-B', '-S', '-m', PREFIX+'actions', str(CASE), 'replay'])]
    dump(CASE/'protocol.json', {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'same_original_H10_lower_dictionary': True, 'fresh_zero_Gram_and_ideal_coordinates': True,
        'proposal_initial_dual': 'Fresh canonical HF product moments',
        'accepting_checker_unchanged': True, 'upper_dependency': str(UPPER),
        'not_a_success_of_the_original_frozen_random_state_rule': True,
        'original_state_failures_remain_additional_cost': True, 'steps': steps})
    started = time.monotonic()
    outcomes = []
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h10_rotated_lower_'+name,
            '--seconds', str(seconds), '--']+command, cwd=ROOT)
        outcomes.append({'stage': name, 'exit_code': result.returncode})
        if result.returncode:
            break
    dump(CASE/'execution.json', {'steps': outcomes, 'elapsed_seconds': time.monotonic()-started,
        'status': 'completed' if len(outcomes) == len(steps) and not outcomes[-1]['exit_code'] else 'stopped_after_failed_step',
        'last_step': name, 'upper_and_integrals_costs_are_additional': True})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'initialize'))
    a = p.parse_args()
    run() if a.action == 'run' else initialize()
