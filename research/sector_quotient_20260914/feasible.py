"""Fixed-bound affine/cone intersection with a sparse orbital preconditioner.

All old compact spans and cross terms remain available. Only the numerical
method changes: seek a positive certificate at a chosen bound directly.
"""
import argparse
from fractions import Fraction as F
import json
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg, splu
from research.sector_quotient_20260914.search import Operator, export
from research.sector_quotient_20260914.eliminated import Quotient
from research.sector_quotient_20260914.boundary import positive
from research.sector_quotient_20260914.budget import OUT
from research.reconstruction_compression_20260914.inputs import dump


def full_preconditioner(op):
    """Drop the reduced-span projections for a sparse normal-equation metric."""
    started = time.monotonic()
    G = (op.free@op.free.T).tocsc()
    for M, V in zip(op.M, op.V):
        n = V.shape[0]
        transposed = (np.arange(n)[None, :]*n+np.arange(n)[:, None]).ravel()
        S = (M+M[:, transposed])*.5
        G = G+S@S.T
    G.eliminate_zeros()
    if G.nnz > 12_000_000:
        raise ValueError('Preconditioner exceeds declared sparse allocation budget')
    regularizer = 1e-10
    factor = splu((G+regularizer*sparse.eye(G.shape[0], format='csc')).tocsc())
    record = {'dimension': G.shape[0], 'normal_nonzeros': G.nnz,
              'factor_nonzeros': factor.L.nnz+factor.U.nnz,
              'diagonal_regularizer': regularizer, 'seconds': time.monotonic()-started,
              'full_reduced_congruence_map_constructed': False,
              'role': 'numerical preconditioner only, no energy acceptance'}
    return factor.solve, record


def run(tag, width_mHa=1.4, seconds=300, max_cg=120, restart=None, normal_cache=None):
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed campaign')
    started = time.monotonic(); prepared = OUT/'prepared'
    folder = OUT/'candidates'/tag; folder.mkdir(parents=True, exist_ok=False)
    meta = json.loads((prepared/'frame.json').read_text()); op = Operator(prepared, meta)
    checks = op.check_adjoint(); quotient = Quotient(op)
    if normal_cache:
        from pathlib import Path
        from scipy.linalg import cho_solve
        from research.reconstruction_compression_20260914.inputs import sha
        cache = Path(normal_cache); preconditioner_record = json.loads((cache/'receipt.json').read_text())
        if (sha(prepared/'frame.json') != preconditioner_record['prepared_frame_sha256'] or
            sha(cache/'cholesky.npy') != preconditioner_record['factor_sha256']):
            raise ValueError('Changed normal preconditioner')
        factor = np.load(cache/'cholesky.npy')
        inverse = lambda y: cho_solve((factor, True), y, check_finite=False)
    else:
        inverse, preconditioner_record = full_preconditioner(op)
    frozen = json.loads((OUT/'frozen_inputs.json').read_text())['h8']
    target = F(frozen['upper_Ha'])-F(str(width_mHa))/1000
    b = float(target)
    def project(y): return y-quotient.U@(quotient.U.T@y)
    def A(Q): return project(op.A(Q))
    def AT(y): return op.AT(project(y))
    rhs = project(op.rhs-op.free[:, 0].toarray().ravel()*b)
    operator = LinearOperator((len(rhs), len(rhs)), matvec=lambda y: A(AT(y)))
    preconditioner = LinearOperator(operator.shape, matvec=lambda y: project(inverse(project(y))))
    Q = [positive(q) for q in op.Q]; z = [q.copy() for q in Q]
    lagrange = np.zeros(len(rhs))
    if restart:
        raw = np.load(restart)
        z = [raw[f'z_{i}'] for i in range(len(z))]
        lagrange = project(raw['lagrange'])
    best = {'residual': float('inf')}; history = []
    details = {'kind': 'full_quotient_fixed_bound_feasibility', 'target_b': str(target),
               'target_width_before_residual_mHa': width_mHa, 'gram_entries': sum(q.size for q in Q),
               'coefficient_rows': len(rhs), 'all_compact_span_cross_terms_retained': True,
               'independent_partner_Gram_blocks': True, 'quartic_number_basis_dimension': 1484,
               'many_body_states_enumerated': 0, 'full_Gram_teacher_used': False, 'additional_MPS_moments_used': False,
               'quotient_conditioning': quotient.record, 'adjoint_check': checks,
               'preconditioner': preconditioner_record, 'setup_seconds': time.monotonic()-started}
    dump(folder/'construction.json', details); print(json.dumps(details), flush=True)
    last_log = time.monotonic()
    def record(iteration, status, count, normal_error):
        x = quotient.recover(Q)
        # Keep the requested b. The difference is visible in the residual.
        x[0] = b
        remaining = op.rhs-op.A(Q)-op.free[:, 0].toarray().ravel()*b
        from scipy.linalg import solve_triangular
        x[quotient.pivots] = solve_triangular(quotient.R, quotient.U.T@remaining)
        residual = (op.A(Q)+op.free@x-op.rhs)/op.scale
        l1 = float(np.sum(abs(residual)))
        rec = {'iteration': iteration, 'seconds': time.monotonic()-started,
               'b': b, 'selected_row_l1': l1, 'affine_residual_l2': float(np.linalg.norm(A(Q)-rhs)),
               'cg_status': status, 'cg_iterations': count, 'cg_actual_relative_residual': normal_error}
        history.append(rec); print(json.dumps(rec), flush=True)
        if l1 < best['residual']:
            best.update(residual=l1, Q=[q.copy() for q in Q], x=x.copy(), record=rec)
        dump(folder/'history.json', history)
        np.savez_compressed(folder/'checkpoint.npz', x=x, lagrange=lagrange,
                            **{f'z_{i}': v for i, v in enumerate(z)}, **{f'Q_{i}': q for i, q in enumerate(Q)})
        return rec
    record(0, 0, 0, 0.)
    solve_start = time.monotonic()
    for iteration in range(10000):
        if time.monotonic()-solve_start >= seconds: break
        Q = [positive(v) for v in z]
        reflected = [2*q-v for q, v in zip(Q, z)]
        residual = A(reflected)-rhs
        counter = [0]
        def count(_): counter[0] += 1
        lagrange, status = cg(operator, residual, x0=lagrange, rtol=1e-7, atol=1e-12,
                              maxiter=max_cg, M=preconditioner, callback=count)
        lagrange = project(lagrange)
        correction = AT(lagrange)
        z = [q-v for q, v in zip(Q, correction)]
        if iteration == 0 or time.monotonic()-last_log >= 10:
            error = float(np.linalg.norm(operator@lagrange-residual)/max(np.linalg.norm(residual), 1e-30))
            rec = record(iteration+1, int(status), counter[0], error); last_log = time.monotonic()
            if rec['selected_row_l1'] < 1e-8: break
    Q = [positive(v) for v in z]
    record('final', int(status), counter[0], float(np.linalg.norm(operator@lagrange-residual)/max(np.linalg.norm(residual), 1e-30)))
    result = export(op, best['Q'], best['x'], folder/'best')
    details.update(best=result, total_seconds=time.monotonic()-started,
                   status='requires_independent_exact_replay', feasibility_proved=False,
                   note='Failure to converge in this finite run is not a cone obstruction.')
    dump(folder/'discovery.json', details)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('tag'); p.add_argument('--width-mHa', type=float, default=1.4)
    p.add_argument('--seconds', type=int, default=300); p.add_argument('--max-cg', type=int, default=120); p.add_argument('--restart')
    p.add_argument('--normal-cache')
    a = p.parse_args(); run(a.tag, a.width_mHa, a.seconds, a.max_cg, a.restart, a.normal_cache)
