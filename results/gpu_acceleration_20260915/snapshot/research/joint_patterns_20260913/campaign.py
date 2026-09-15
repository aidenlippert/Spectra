"""Serial, process-isolated measured campaign with an outer wall watchdog."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CASES = [('h4', 6, 0, 'joint'), ('h4', 6, 6, 'joint'),
    ('h6', 10, 0, 'joint'), ('h6', 10, 10, 'separate'),
    ('h6', 10, 2, 'joint'), ('h6', 10, 4, 'joint'), ('h6', 10, 10, 'joint')]


def run(out, budget=90.):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    prior = ROOT/'results/molecular_collective_20260913/campaign'
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    results = []
    for name, rank, count, coupling in CASES:
        label = f'{name}_{coupling}_{count}'; directory = out/label
        fixture = prior/name/'fixture.json'; tail = prior/name/f'rank_{rank}'/'tail.json'
        cmd = [sys.executable, '-m', 'research.joint_patterns_20260913.discovery',
            '--fixture', str(fixture), '--tail', str(tail), '--out', str(directory),
            '--count', str(count), '--coupling', coupling, '--budget', str(budget)]
        before = time.monotonic(); failure = None
        with (out/f'{label}.log').open('w') as log:
            try:
                proc = subprocess.run(cmd, cwd=ROOT, env=environment, stdout=log, stderr=subprocess.STDOUT, timeout=budget)
                if proc.returncode:
                    failure = f'exit_{proc.returncode}'
            except subprocess.TimeoutExpired:
                failure = 'end_to_end_wall_timeout'
        wall = time.monotonic()-before
        accepted = failure is None and wall <= budget and (directory/'receipt.json').exists()
        row = {'case': label, 'system': name, 'constraint_patterns': count, 'coupling': coupling,
            'wall_seconds': wall, 'budget_seconds': budget, 'within_budget_accepted': accepted, 'failure': failure,
            'fixture_path': str(fixture.relative_to(ROOT)), 'tail_path': str(tail.relative_to(ROOT)),
            'fixture_file_sha256': hashlib.sha256(fixture.read_bytes()).hexdigest(),
            'tail_file_sha256': hashlib.sha256(tail.read_bytes()).hexdigest()}
        if accepted:
            receipt = json.loads((directory/'receipt.json').read_text())
            row.update({'original_lower_Ha': receipt['accepted']['original_lower_Ha'],
                'original_lower_float_Ha': receipt['accepted']['original_lower_float_Ha'],
                'compact_certificate_bytes': receipt['compact_certificate_bytes'], 'solver_status': receipt['status']})
        results.append(row)
        (out/'summary.json').write_text(json.dumps({'cases': results, 'wall_seconds': time.monotonic()-start}, indent=2)+'\n')
        print(json.dumps(row), flush=True)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--budget', type=float, default=90.)
    args = parser.parse_args(); run(args.out, args.budget)
