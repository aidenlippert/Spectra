"""Synchronized FP64 PSD projection benchmarks on the real block dimensions."""
import argparse
import json
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_limits


def cpu_project(blocks, driver='evr'):
    result = []
    for q in blocks:
        if not np.isfinite(q).all():
            raise ValueError('Nonfinite PSD input')
        kwargs = dict(driver='evr', subset_by_value=(0., np.inf)) if driver == 'positive' else dict(driver=driver)
        e, v = linalg.eigh((q + q.T) * .5, check_finite=False, **kwargs)
        result.append((v * np.maximum(e, 0)) @ v.T)
    return result


class GPUProject:
    def __init__(self, shapes, batched=True, mixed=False):
        import cupy as cp
        self.cp = cp
        self.shapes = list(shapes)
        groups = defaultdict(list)
        for i, shape in enumerate(shapes):
            groups[shape].append(i)
        self.groups = list(groups.values()) if batched else [[i] for i in range(len(shapes))]
        if mixed:
            self.groups = [part for ids in self.groups for part in
                           ([[i] for i in ids] if shapes[ids[0]][0] > 64 else [ids])]

    def device(self, blocks):
        cp = self.cp
        result = [None] * len(blocks)
        for ids in self.groups:
            if len(ids) == 1:
                q = blocks[ids[0]]
                e, v = cp.linalg.eigh((q + q.T) * .5)
                result[ids[0]] = (v * cp.maximum(e, 0)) @ v.T
                continue
            q = cp.stack([blocks[i] for i in ids])
            q = (q + q.swapaxes(-1, -2)) * .5
            e, v = cp.linalg.eigh(q)
            projected = (v * cp.maximum(e, 0)[:, None, :]) @ v.swapaxes(-1, -2)
            for k, i in enumerate(ids):
                result[i] = projected[k]
        return result

    def __call__(self, blocks):
        cp = self.cp
        if [q.shape for q in blocks] != self.shapes:
            raise ValueError('Changed projection block shapes')
        if not all(np.isfinite(q).all() for q in blocks):
            raise ValueError('Nonfinite PSD input')
        offsets = np.cumsum([0] + [q.size for q in blocks])
        packed = cp.asarray(np.concatenate([q.ravel() for q in blocks]))
        dev = [packed[offsets[i]:offsets[i+1]].reshape(q.shape) for i, q in enumerate(blocks)]
        import cupyx
        with cupyx.errstate(linalg='raise'):
            answer = self.device(dev)
        packed_out = cp.asnumpy(cp.concatenate([q.ravel() for q in answer]))
        if not np.isfinite(packed_out).all():
            raise ValueError('Nonfinite GPU projection')
        return [packed_out[offsets[i]:offsets[i+1]].reshape(q.shape) for i, q in enumerate(blocks)]


class HybridProject:
    """Leave small blocks on CPU; send only larger eigensolves to CUDA."""
    def __init__(self, shapes, cutoff=128):
        self.shapes = list(shapes)
        self.small = [i for i, shape in enumerate(shapes) if shape[0] < cutoff]
        self.large = [i for i, shape in enumerate(shapes) if shape[0] >= cutoff]
        self.gpu = GPUProject([shapes[i] for i in self.large], batched=False) if self.large else None

    def __call__(self, blocks):
        if [q.shape for q in blocks] != self.shapes:
            raise ValueError('Changed projection block shapes')
        result = [None] * len(blocks)
        for i, q in zip(self.small, cpu_project([blocks[i] for i in self.small], 'evd')):
            result[i] = q
        if self.gpu:
            for i, q in zip(self.large, self.gpu([blocks[i] for i in self.large])):
                result[i] = q
        return result


def bench(case, output, checkpoint=None):
    import cupy as cp
    z = np.load(case / 'prepared/bases.npz')
    sizes = [z[f'V_{i}'].shape[1] for i in range(len(z.files))]
    rng = np.random.default_rng(51)
    blocks = [rng.normal(size=(n, n)) for n in sizes]
    blocks = [(q + q.T) * .5 for q in blocks]
    if checkpoint:
        from scipy import sparse
        state = np.load(checkpoint)
        blocks = [state[f'Q_{i}'] + .03 * (sparse.load_npz(case/f'prepared/map_{i}.npz').T @ state['y']).reshape(n,n)
                  for i,n in enumerate(sizes)]
    with threadpool_limits(limits=1):
        reference = cpu_project(blocks)
    results = []

    def trial(name, fn, sync=lambda: None, threads=1):
        with threadpool_limits(limits=threads):
            got = fn()
            sync()
            error = max(np.linalg.norm(a - b) / max(1., np.linalg.norm(b))
                        for a, b in zip(got, reference))
            if error > 1e-10:
                raise ValueError((name, error))
            times = []
            for _ in range(7):
                sync()
                start = time.perf_counter()
                for _ in range(5):
                    fn()
                sync()
                times.append((time.perf_counter() - start) / 5)
        row = dict(name=name, threads=threads, median_seconds=float(np.median(times)),
                   min_seconds=min(times), max_seconds=max(times), relative_error=error)
        results.append(row)
        print(json.dumps(row), flush=True)

    # Higher BLAS thread counts are tested in separate processes: the remote
    # OpenBLAS build crashed while changing counts in-process in the first probe.
    for threads in (1,):
        for driver in ('evr', 'evd', 'positive'):
            trial(f'cpu_{driver}', lambda d=driver: cpu_project(blocks, d), threads=threads)
    for workers in (2, 4, 8):
        with ThreadPoolExecutor(workers) as pool:
            trial(f'cpu_evd_independent_blocks_{workers}_workers',
                  lambda: list(pool.map(lambda q: cpu_project([q], 'evd')[0], blocks)))
    for cutoff in (64, 128, 200):
        hybrid = HybridProject([q.shape for q in blocks], cutoff)
        trial(f'hybrid_cutoff_{cutoff}', lambda: hybrid(blocks), cp.cuda.Stream.null.synchronize)
    for label, batched, mixed in [('individual', False, False), ('batched', True, False), ('mixed', True, True)]:
        gpu = GPUProject([q.shape for q in blocks], batched, mixed)
        trial(f'gpu_{label}_including_transfers',
              lambda: gpu(blocks), cp.cuda.Stream.null.synchronize)
        dev = [cp.asarray(q) for q in blocks]
        # Correctness is checked above. Only output copies are excluded in this timing.
        gpu.device(dev)
        cp.cuda.Stream.null.synchronize()
        times = []
        for _ in range(7):
            start = time.perf_counter()
            for _ in range(5):
                gpu.device(dev)
            cp.cuda.Stream.null.synchronize()
            times.append((time.perf_counter() - start) / 5)
        row = dict(name=f'gpu_{label}_resident',
                   median_seconds=float(np.median(times)), min_seconds=min(times), max_seconds=max(times))
        results.append(row)
        print(json.dumps(row), flush=True)
    output.write_text(json.dumps(dict(case=str(case), sizes=sizes, numpy=np.__version__,
                                     cupy=cp.__version__, gpu=cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
                                     results=results), indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--checkpoint', type=Path)
    a = p.parse_args()
    bench(a.case, a.output, a.checkpoint)
