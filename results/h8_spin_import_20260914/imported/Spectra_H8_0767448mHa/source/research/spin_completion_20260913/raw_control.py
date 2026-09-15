"""One serial bounded full-family control without numerical whitening."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time

from research.spin_completion_20260913.campaign import ROOT, save


def arm(out, start, budget):
    from research.spin_completion_20260913.discovery import Model
    out.mkdir(parents=True, exist_ok=False)
    prior = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((prior/'fixture.json').read_text()); tail = json.loads((prior/'rank_10/tail.json').read_text())
    model = Model(data, tail)
    summary = {'arm': 'full_raw', 'budget_seconds': budget, 'model_construction': model.construction,
        'complete': False, 'trials': [], 'best_case': None}
    save(out/'summary.json', summary)
    before = time.monotonic()
    try:
        receipt, _ = model.solve(model.full_span(), out/'full_frame', budget-(time.monotonic()-start), solver_cap=150., whiten=False)
        row = {'case': 'full_frame', 'accepted': True, 'directions': 492,
            'lower_Ha': receipt['accepted']['original_lower_Ha'], 'trial_wall_seconds': time.monotonic()-before,
            'complete_at_seconds': time.monotonic()-start, 'certificate_bytes': receipt['certificate_bytes']}
        if row['complete_at_seconds'] <= budget:
            summary['best_case'] = 'full_frame'
    except Exception as error:
        row = {'case': 'full_frame', 'accepted': False, 'directions': 492, 'error': repr(error),
            'trial_wall_seconds': time.monotonic()-before, 'complete_at_seconds': time.monotonic()-start}
    summary.update({'complete': True, 'trials': [row], 'elapsed_seconds': time.monotonic()-start})
    save(out/'summary.json', summary); print(json.dumps(row), flush=True)


def run(out):
    if out.exists() or out.with_suffix('.log').exists() or out.with_name(out.name+'_watchdog.json').exists():
        raise FileExistsError('Refusing to overwrite a full-control run or its sidecars')
    start = time.monotonic(); budget = 240.; failure = None
    environment = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    with out.with_suffix('.log').open('w') as log:
        try:
            proc = subprocess.run([sys.executable, '-m', 'research.spin_completion_20260913.raw_control', '--out', str(out),
                '--started', str(start), '--budget', str(budget)], env=environment, cwd=ROOT,
                stdout=log, stderr=subprocess.STDOUT, timeout=budget)
            if proc.returncode:
                failure = f'exit_{proc.returncode}'
        except subprocess.TimeoutExpired:
            failure = 'end_to_end_wall_timeout'
    row = {'arm': 'full_raw', 'budget_seconds': budget, 'wall_seconds': time.monotonic()-start, 'failure': failure}
    if (out/'summary.json').exists():
        summary = json.loads((out/'summary.json').read_text())
        row.update({'complete': summary['complete'], 'best_case': summary['best_case'],
            'completed_valid_trials': sum(r['accepted'] and r['complete_at_seconds'] <= budget for r in summary['trials'])})
    save(out.with_name(out.name+'_watchdog.json'), row); print(json.dumps(row), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--started', type=float); parser.add_argument('--budget', type=float, default=240.)
    args = parser.parse_args()
    if args.started is None:
        run(args.out)
    else:
        arm(args.out, args.started, args.budget)
