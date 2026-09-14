"""Measured, bounded float64 search throughput on frozen physical matrices."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.linalg import eigvalsh


def evaluate_groups(xp, groups, delta):
    minimum = xp.full(len(delta), xp.inf, dtype=xp.float64)
    for n, base, derivative in groups:
        matrices = (delta @ derivative).reshape(len(delta), len(base), n, n) + base
        values = matrices[..., 0, 0] if n == 1 else xp.linalg.eigvalsh(matrices.reshape(-1, n, n))[:, 0].reshape(len(delta), len(base))
        minimum = xp.minimum(minimum, values.min(axis=1))
    return minimum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--evaluations', type=int, default=256)
    ap.add_argument('--batch', type=int, default=128)
    args = ap.parse_args()
    if not 1 <= args.evaluations <= 1024 or not 1 <= args.batch <= 256:
        raise ValueError('Bounded to 1024 evaluations and batch 256')
    start = time.perf_counter()
    data = np.load(args.data, allow_pickle=False)
    indices = sorted(int(k[1:]) for k in data.files if k.startswith('A'))
    by_size = {}
    scalar = []
    for i in indices:
        base = data[f'A{i}'] + (float(data['alpha']) + float(data['beta'])) * data[f'P{i}'] + float(data['beta']) * float(data['ratio']) * data[f'Q{i}']
        derivative = data[f'D{i}']
        scalar.append((base, derivative))
        by_size.setdefault(len(base), []).append((base, derivative))
    groups = []
    for n, items in sorted(by_size.items()):
        base = np.stack([a for a, _ in items])
        derivative = np.stack([d for _, d in items]).transpose(1, 0, 2, 3).reshape(6, -1)
        groups.append((n, base, derivative))
    delta = np.random.default_rng(20260912).normal(0, .002, (args.evaluations, 6))
    delta[0] = 0
    load_seconds = time.perf_counter() - start
    start = time.perf_counter()
    reference = np.array([min(float(eigvalsh(a + np.einsum('i,ijk->jk', row, d), subset_by_index=[0, 0], check_finite=False)[0]) for a, d in scalar) for row in delta])
    scalar_seconds = time.perf_counter() - start
    start = time.perf_counter()
    cpu = np.concatenate([evaluate_groups(np, groups, delta[j:j+args.batch]) for j in range(0, len(delta), args.batch)])
    cpu_seconds = time.perf_counter() - start
    if not np.allclose(cpu, reference, atol=1e-10, rtol=1e-10):
        raise ValueError('Grouped CPU disagrees with independent scalar assembly')
    import cupy as cp
    import cupyx
    cp.cuda.Device(0).use()
    sync = cp.cuda.Stream.null.synchronize
    start = time.perf_counter()
    device_groups = [(n, cp.asarray(a), cp.asarray(d)) for n, a, d in groups]
    device_delta = cp.asarray(delta)
    sync()
    transfer_seconds = time.perf_counter() - start
    # Report a full-shape warmup separately, including compilation and library initialization.
    start = time.perf_counter()
    with cupyx.errstate(linalg='raise'):
        evaluate_groups(cp, device_groups, device_delta[:args.batch])
        sync()
    warmup_seconds = time.perf_counter() - start
    timings = []
    for repeat in range(3):
        start = time.perf_counter()
        with cupyx.errstate(linalg='raise'):
            results = [evaluate_groups(cp, device_groups, device_delta[j:j+args.batch]) for j in range(0, len(delta), args.batch)]
            got = cp.asnumpy(cp.concatenate(results))
            sync()
        timings.append(time.perf_counter() - start)
        if not np.allclose(got, reference, atol=1e-9, rtol=1e-10):
            raise ValueError('GPU differs from scalar reference')
    median = float(np.median(timings))
    props = cp.cuda.runtime.getDeviceProperties(0)
    name = props['name'].decode() if isinstance(props['name'], bytes) else props['name']
    receipt = {
        'scope': 'Numerical physical-matrix batch search, not an exact certificate',
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'data_sha256': hashlib.sha256(args.data.read_bytes()).hexdigest(),
        'evaluations': args.evaluations, 'batch': args.batch, 'sectors': len(indices),
        'shape_groups': len(groups), 'load_and_group_seconds': load_seconds,
        'cpu_scalar_seconds': scalar_seconds, 'cpu_grouped_seconds': cpu_seconds,
        'gpu_transfer_seconds': transfer_seconds, 'gpu_warmup_seconds': warmup_seconds,
        'gpu_seconds_repeats_including_result_copy': timings,
        'gpu_median_seconds': median,
        'gpu_cold_seconds_excluding_matrix_generation': load_seconds + transfer_seconds + warmup_seconds + timings[0],
        'resident_speedup_vs_faster_measured_cpu': min(scalar_seconds, cpu_seconds) / median,
        'gpu_max_abs_error': float(np.max(np.abs(got-reference))),
        'cpu_grouped_max_abs_error': float(np.max(np.abs(cpu-reference))),
        'validated': True, 'gpu_name': name,
        'resident_bytes': sum(a.nbytes + d.nbytes for _, a, d in groups) + delta.nbytes,
        'max_matrix_batch_bytes': max(args.batch * len(a) * n * n * 8 for n, a, _ in groups),
        'versions': {'numpy': np.__version__, 'cupy': cp.__version__},
        'cpu_values': reference.tolist(), 'gpu_values': got.tolist(),
    }
    args.out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: v for k, v in receipt.items() if not k.endswith('_values')}, indent=2))


if __name__ == '__main__':
    main()
