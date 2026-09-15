"""Measured library candidates on retained Spectra arrays; no accepting claims."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import traceback


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2)+'\n')


def sparse_bench(case, output):
    import numpy as np
    from scipy import sparse
    from scipy.sparse.linalg import splu
    from sksparse.cholmod import cholesky
    from threadpoolctl import threadpool_limits
    import cupy as cp
    from cupyx.scipy.sparse import csr_matrix
    from nvmath.sparse.advanced import DirectSolver, DirectSolverMatrixType
    start = time.perf_counter()
    path = case/'prepared/normal.npz'
    G = sparse.load_npz(path).tocsc()
    A = G+sparse.eye(G.shape[0], format='csc')*1e-11
    rng = np.random.default_rng(109)
    # Reproducible, consistent RHS; the coefficient matrix is the actual normal
    # system used by the solver. This is not a complete optimization benchmark.
    b = np.asarray(A@rng.normal(size=A.shape[0]))
    record = {'case': str(case), 'normal_file_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'shape': A.shape, 'nnz': A.nnz, 'diagonal_shift': 1e-11, 'precision': 'float64',
        'blas_threads': 1, 'RHS': 'Actual shifted normal matrix times deterministic seed-109 vector',
        'input_seconds': time.perf_counter()-start, 'results': [], 'full_solver_speedup_measured': False}
    save(output, record)

    def measure(name, setup):
        row = {'backend': name}
        try:
            t = time.perf_counter()
            solve, close = setup()
            cp.cuda.Stream.null.synchronize()
            row['setup_seconds'] = time.perf_counter()-t
            got = solve(b)
            times = []
            for _ in range(5):
                t = time.perf_counter()
                for _ in range(20):
                    got = solve(b)
                cp.cuda.Stream.null.synchronize()
                times.append((time.perf_counter()-t)/20)
            r = A@got-b
            row.update(median_solve_seconds=float(np.median(times)), min_solve_seconds=min(times),
                max_solve_seconds=max(times), relative_RHS_residual=float(np.linalg.norm(r)/np.linalg.norm(b)),
                backward_error=float(np.linalg.norm(r)/(np.linalg.norm(A.data)*np.linalg.norm(got)+np.linalg.norm(b))),
                status='passed' if np.isfinite(got).all() and np.linalg.norm(r) <= 1e-8*np.linalg.norm(b) else 'accuracy_gate_failed')
            close()
        except Exception as exc:
            row.update(status='failed', error=repr(exc), traceback=traceback.format_exc())
        record['results'].append(row)
        save(output, record)
        print(json.dumps(row), flush=True)

    def cpu_lu():
        factor = splu(A)
        return factor.solve, lambda: None

    def cholmod():
        factor = cholesky(A)
        return lambda rhs: np.asarray(factor(rhs)), lambda: None

    def cudss():
        a = csr_matrix(A)
        rhs = cp.asarray(b)
        solver = DirectSolver(a, rhs, options={'sparse_system_type': DirectSolverMatrixType.SPD})
        solver.plan()
        solver.factorize()
        def solve(value):
            rhs.set(value)
            return cp.asnumpy(solver.solve())
        return solve, solver.free

    with threadpool_limits(limits=1):
        measure('SciPy_SuperLU', cpu_lu)
        measure('SuiteSparse_CHOLMOD', cholmod)
        measure('NVIDIA_cuDSS_with_RHS_and_output_transfers', cudss)
    record['complete_seconds'] = time.perf_counter()-start
    save(output, record)


def tensor_bench(case, output):
    import numpy as np
    import cupy as cp
    from cuquantum.tensornet import Network
    from research.reconstruction_compression_20260914.moments import Proposal, transfer
    from research.correlated_pair_20260913.mps_exact import I
    data = json.loads((case/'fixture.json').read_text())
    state = json.loads((case/'mps/state.json').read_text())
    start = time.perf_counter()
    oracle = Proposal(data, state)
    site = max(range(oracle.m), key=lambda i: oracle.dims[i]*oracle.dims[i+1])
    E = np.ones((1, 1))
    for i in range(site):
        E = transfer(E, oracle.edges[i], np.asarray(I), oracle.dims[i+1])
    edges = oracle.edges[site]
    T = np.zeros((oracle.dims[site], 2, oracle.dims[site+1]))
    for a, s, b, v in edges:
        T[int(a), int(s), int(b)] = v
    op = np.array([[1., .25], [.25, -1.]])
    operands = [E, T, T, op]
    expression = 'ab,asr,btq,st->rq'
    reference = np.einsum(expression, *operands, optimize=True)
    result = {'case': str(case), 'site': site, 'tensor_shapes': [a.shape for a in operands],
        'source_state_sha256': hashlib.sha256((case/'mps/state.json').read_bytes()).hexdigest(),
        'precision': 'float64', 'expression': expression, 'input_and_JIT_seconds': time.perf_counter()-start,
        'scope': 'One actual MPS transfer, not complete DMRG discovery', 'results': []}
    save(output, result)
    def measure(name, fn, gpu=False):
        got = fn()
        if isinstance(got, cp.ndarray):
            got = cp.asnumpy(got)
        error = float(np.linalg.norm(got-reference)/max(np.linalg.norm(reference),1e-30))
        if error > 1e-10:
            raise ValueError((name,'contraction disagreement',error))
        times = []
        for _ in range(5):
            if gpu:
                cp.cuda.Stream.null.synchronize()
            t = time.perf_counter()
            for _ in range(50):
                fn()
            if gpu:
                cp.cuda.Stream.null.synchronize()
            times.append((time.perf_counter()-t)/50)
        row = {'backend': name, 'median_seconds': float(np.median(times)), 'relative_error': error,
            'min_seconds': min(times), 'max_seconds': max(times)}
        result['results'].append(row)
        save(output, result)
        print(json.dumps(row),flush=True)
    measure('existing_NumPy_numba_sparse_transfer', lambda: transfer(E, edges, op.ravel(), oracle.dims[site+1]))
    measure('NumPy_dense_einsum', lambda: np.einsum(expression,*operands,optimize=True))
    t = time.perf_counter()
    dev = [cp.asarray(a) for a in operands]
    cp.cuda.Stream.null.synchronize()
    result['initial_input_transfer_seconds'] = time.perf_counter()-t
    measure('CuPy_dense_einsum_resident',lambda:cp.einsum(expression,*dev,optimize=True),True)
    t = time.perf_counter()
    with Network(expression,*dev) as network:
        network.contract_path()
        cp.cuda.Stream.null.synchronize()
        result['cuTensorNet_plan_seconds'] = time.perf_counter()-t
        measure('cuTensorNet_resident_cached_path',lambda:network.contract(),True)
        def transferred():
            arrays = [cp.asarray(a) for a in operands]
            network.reset_operands(*arrays)
            return cp.asnumpy(network.contract())
        measure('cuTensorNet_with_input_and_output_transfers',transferred,True)
    result['complete_seconds'] = time.perf_counter()-start
    save(output,result)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('mode',choices=('sparse','tensor'))
    p.add_argument('case',type=Path)
    p.add_argument('output',type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise ValueError('Preserve existing benchmark results')
    {'sparse':sparse_bench,'tensor':tensor_bench}[a.mode](a.case.resolve(),a.output.resolve())
