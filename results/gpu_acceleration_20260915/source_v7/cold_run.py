"""Measure the selected acceleration in a fresh, isolated molecular calculation."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def run(case, state_python=None, fast_prepare=False):
    case = case.resolve()
    case.mkdir(parents=True, exist_ok=False)
    spec = {'geometry': [['H', [0.07*i*i, 0.03*(i % 3), 1.4*i]] for i in range(6)],
            'basis': 'sto-3g', 'role': 'fresh asymmetric H6 acceleration measurement'}
    (case/'specification.json').write_text(json.dumps(spec, indent=2)+'\n')
    prefix = 'research.transfer_solver_20260915.'
    py = sys.executable
    mod = lambda name, args, exact=False: [py, '-B']+(['-S'] if exact else [])+['-m', name]+args
    c = str(case)
    steps = [
        ('generate', mod(prefix+'generate', [c]), 90),
        ('state', mod(prefix+'actions', [c, 'state']), 600),
        ('upper', mod(prefix+'actions', [c, 'upper'], True), 180),
        ('nonsinglet', mod(prefix+'actions', [c, 'nonsinglet']), 240),
        ('prepare', mod(prefix+'prepare', [c]), 600),
        ('optimize', mod('research.gpu_acceleration_20260915.solve',
                         [c, 'hybrid_adaptive', '--seconds', '240', '--mu', '2', '--backend', 'hybrid', '--adaptive']), 300),
        ('replay', mod('research.gpu_acceleration_20260915.replay', [c, 'hybrid_adaptive'], True), 300)]
    if state_python:
        # Preserve the virtual-environment entry point rather than dereferencing
        # its interpreter symlink and losing its installed dependencies.
        steps[1][1][0] = str(state_python.absolute())
    if fast_prepare:
        steps[4] = ('prepare', mod('research.gpu_acceleration_20260915.prepare', [c]), 600)
    source = Path(__file__).resolve().parents[2]
    hashes = {str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
              for name in ('transfer_solver_20260915','gpu_acceleration_20260915')
              for p in (source/'research'/name).glob('*.py')}
    protocol = dict(created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    specification=spec, source_hashes=hashes,
                    backend='A100 hybrid FP64; CPU blocks smaller than 128',
                    refinement_rule='mu 2 -> .03 when predicted width <10 mHa and primal residual <1e-4, checked every100 iterations',
                    numerical_stop='predicted width <1.3 mHa and primal residual <5e-8',
                    exact_acceptance='unchanged rational upper and complete-spin lower checkers; width <=1.6mHa',
                    historical_states_maps_checkpoints_or_proofs_loaded=False,
                    existing_fast_exact_projector_used=fast_prepare,
                    dependencies_installation_and_machine_boot_excluded=True,
                    steps=[dict(name=n, command=cmd, timeout=t) for n,cmd,t in steps])
    (case/'protocol.json').write_text(json.dumps(protocol, indent=2)+'\n')
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMBA_CACHE_DIR=str(case/'numba_cache'), PYTHONDONTWRITEBYTECODE='1',
               PYTHONPATH=str(source))
    results = []
    started = time.monotonic()
    for name, cmd, timeout in steps:
        for p, expected in hashes.items():
            if hashlib.sha256((source/p).read_bytes()).hexdigest() != expected:
                raise ValueError('Frozen campaign source changed')
        begin = time.monotonic()
        with (case/f'{name}.log').open('x') as log:
            try:
                result = subprocess.run(cmd, cwd=source, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
                status = result.returncode
            except subprocess.TimeoutExpired:
                status = 'timeout'
        row = dict(step=name, seconds=time.monotonic()-begin, exit_code=status)
        results.append(row)
        (case/'execution.json').write_text(json.dumps(dict(steps=results, total_seconds=time.monotonic()-started,
                                                    completed=len(results)==len(steps) and status==0), indent=2)+'\n')
        print(json.dumps(row), flush=True)
        if status != 0:
            raise SystemExit(1)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('--state-python', type=Path)
    p.add_argument('--fast-prepare', action='store_true')
    args = p.parse_args()
    run(args.case, args.state_python, args.fast_prepare)
