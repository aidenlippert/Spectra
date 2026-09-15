"""Named caller for optional sparse-library proposal backends; exact replay separate."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def run(args):
    from research.gpu_acceleration_20260915 import solve as base
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    base.Quotient = SparseQuotient
    calls = [0]
    if args.normal in ('cholmod','cudss'):
        from scipy import sparse
        original = base.splu
        G = sparse.load_npz(args.case/'prepared/normal.npz').tocsc()
        def selected_factor(matrix):
            calls[0] += 1
            if calls[0] > 1:
                return original(matrix)
            expected = G+sparse.eye(G.shape[0],format='csc')*1e-11
            difference = matrix-expected
            if difference.nnz and max(abs(difference.data)) > 1e-15:
                raise ValueError('Unexpected normal-system factorization caller')
            if args.normal == 'cholmod':
                from sksparse.cholmod import cholesky
                factor = cholesky(matrix)
                inverse = lambda rhs: factor(rhs)
            else:
                import cupy as cp
                from cupyx.scipy.sparse import csr_matrix
                from nvmath.sparse.advanced import DirectSolver, DirectSolverMatrixType
                buffer = cp.zeros(matrix.shape[0], dtype=cp.float64)
                factor = DirectSolver(csr_matrix(matrix),buffer,
                    options={'sparse_system_type':DirectSolverMatrixType.SPD})
                factor.plan()
                factor.factorize()
                def inverse(rhs):
                    buffer.set(rhs)
                    return cp.asnumpy(factor.solve())
            class Adapter:
                def solve(self,rhs):
                    return inverse(rhs)
            return Adapter()
        base.splu = selected_factor
    start = time.monotonic()
    base.run(args.case,args.tag,args.seconds,2.,None,args.backend,args.iterations,100,None,args.iterations is None)
    if args.normal in ('cholmod','cudss') and calls[0] != 2:
        raise ValueError('Sparse factorization call pattern changed; audit adapter')
    value = {'sparse_quotient': True, 'normal_factorization': args.normal,
        'seconds': time.monotonic()-start, 'base_solver_source_sha256': hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
        'all_original_solver_residual_and_stopping_gates_preserved': True,
        'fixed_iterations': args.iterations, 'adaptive_schedule': args.iterations is None,
        'exact_acceptance_required': True}
    with (args.case/args.tag/'sparse_backend.json').open('x') as stream:
        json.dump(value,stream,indent=2)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('case',type=Path)
    p.add_argument('tag')
    p.add_argument('--seconds',type=float,required=True)
    p.add_argument('--backend',choices=('cpu_evd','hybrid'),default='hybrid')
    p.add_argument('--normal',choices=('superlu','cholmod','cudss'),default='superlu')
    p.add_argument('--iterations',type=int)
    a=p.parse_args()
    a.case=a.case.resolve()
    run(a)
