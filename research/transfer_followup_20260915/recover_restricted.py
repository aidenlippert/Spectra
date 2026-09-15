"""Export the preserved restricted-family checkpoint after its diagnostic failed.

The numerical reconstruction guard is unchanged. Acceptance recomputes the
complete operator residual with the original integer/rational checker.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
from fractions import Fraction as F
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, PREFIX

CASE = OUT/'controls/h8_restricted64'


def export_saved():
    import numpy as np
    from research.transfer_solver_20260915.solve import Operator
    from research.sector_quotient_20260914.eliminated import Quotient
    from research.sector_quotient_20260914.search import export
    started = time.monotonic()
    checkpoint = CASE/'stage1/checkpoint.npz'
    op = Operator(CASE)
    z = np.load(checkpoint)
    Q = [z[f'Q_{i}'] for i in range(len(op.Q))]
    if any(a.shape != b.shape for a, b in zip(Q, op.Q)):
        raise ValueError('Checkpoint dimensions changed')
    quotient = Quotient(op)
    result = export(op, Q, quotient.recover(Q), CASE/'recovered/export')
    dump(CASE/'recovered/construction.json', {
        'checkpoint_sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        'original_failure_preserved': True, 'new_optimization_iterations': 0,
        'original_numerical_guard_unchanged': True,
        'candidate_only_requires_original_exact_checker': True,
        'export': result, 'seconds': time.monotonic()-started})


def verify():
    import sys
    from research.collective_completion_20260914.spin_screen import check
    read = lambda name: json.loads((CASE/name).read_text())
    start = time.monotonic()
    receipt = check(read('fixture.json'), read('recovered/export/certificate.json'), read('nonsinglet.json'))
    U, L = F(read('upper.json')['upper_Ha']), F(receipt['lower'])
    if L > U:
        raise ValueError('Inconsistent exact endpoints')
    forbidden = [n for n in ('numpy', 'scipy', 'numba', 'quimb', 'pyscf', 'cvxpy') if n in sys.modules]
    if forbidden:
        raise ValueError(('Numerical accepting imports', forbidden))
    dump(CASE/'recovered/lower.json', receipt)
    dump(CASE/'recovered/interval.json', {
        'upper_Ha': str(U), 'lower_Ha': str(L), 'width_Ha': str(U-L),
        'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625),
        'valid_on': receipt['valid_on'], 'no_family_ceiling_claimed': True,
        'accepting_numerical_imports': forbidden, 'seconds': time.monotonic()-start})


def run():
    folder = CASE/'recovered'
    folder.mkdir(exist_ok=False)
    steps = [('export', 90, [NUM, '-B', '-m', 'research.transfer_followup_20260915.recover_restricted', 'export']),
             ('verify', 600, [STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.recover_restricted', 'verify'])]
    dump(folder/'protocol.json', {'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'purpose': 'Recover a concrete exact interval from the last preserved iterate without weakening the numerical guard',
        'new_optimization': False, 'steps': steps})
    outcomes = []
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h8_restricted64_recover_'+name,
            '--seconds', str(seconds), '--']+command, cwd=ROOT)
        outcomes.append({'stage': name, 'exit_code': result.returncode})
        if result.returncode:
            break
    dump(folder/'execution.json', outcomes)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'export', 'verify'))
    a = p.parse_args()
    {'run': run, 'export': export_saved, 'verify': verify}[a.action]()
