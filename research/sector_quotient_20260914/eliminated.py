"""Eliminate number/spin multiplier coordinates before a matrix-free SDP.

This is numerical conditioning of the proposal problem, not deletion of the
contraction-free constraints. Every final factor is checked in the full algebra.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import LinearOperator, cg
from research.sector_quotient_20260914.search import Operator, export
from research.sector_quotient_20260914.boundary import positive
from research.sector_quotient_20260914.budget import ROOT, OUT
from research.reconstruction_compression_20260914.inputs import dump


class Quotient:
    def __init__(self, op):
        started = time.monotonic()
        self.op = op
        F = op.free.toarray()
        # Invariant ideals contain duplicate coordinates. Rank-revealing QR
        # supplies a numerical basis; exact replay accounts for all errors.
        U, R, pivots = linalg.qr(F[:, 1:], mode='economic', pivoting=True)
        diagonal = abs(np.diag(R)); rank = int(np.sum(diagonal > 1e-11*max(diagonal.max(), 1.)))
        self.U = U[:, :rank].copy(); self.R = R[:rank, :rank].copy(); self.pivots = pivots[:rank]+1
        v = F[:, 0]-self.U@(self.U.T@F[:, 0])
        self.y0 = v/np.dot(v, v)
        self.v = v/np.linalg.norm(v)
        self.C = op.AT(self.y0)
        self.offset = float(op.rhs@self.y0)
        self.rhs = self.project(op.rhs)
        self.record = {'ideal_coordinate_count': F.shape[1]-1, 'numerical_ideal_rank': rank,
                       'QR_rank_threshold': 1e-11, 'dependent_coordinates_removed': F.shape[1]-1-rank,
                       'ideal_projector_max_defect': float(np.max(abs(F[:, 1:]-self.U@(self.U.T@F[:, 1:])))),
                       'dual_normalization_error': float(max(abs(F.T@self.y0-np.r_[1., np.zeros(F.shape[1]-1)]))),
                       'seconds': time.monotonic()-started}
        if self.record['ideal_projector_max_defect'] > 1e-8 or self.record['dual_normalization_error'] > 1e-8:
            raise ValueError('Ill-conditioned numerical ideal elimination')

    def project(self, y):
        return y-self.U@(self.U.T@y)-self.v*np.dot(self.v, y)

    def A(self, Q):
        return self.project(self.op.A(Q))

    def AT(self, y):
        return self.op.AT(self.project(y))

    def recover(self, Q):
        residual = self.op.rhs-self.op.A(Q)
        b = float(self.y0@residual)
        values = np.zeros(self.op.free.shape[1]); values[0] = b
        remaining = residual-self.op.free[:, 0].toarray().ravel()*b
        z = self.U.T@remaining
        values[self.pivots] = linalg.solve_triangular(self.R, z)
        return values


def run(tag, seconds=300, max_cg=100, mu=.1, restart=None, normal_cache=None, fixed_mu=False, prepared_path=None):
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed campaign')
    started = time.monotonic(); folder = OUT/'candidates'/tag; folder.mkdir(parents=True, exist_ok=False)
    prepared = Path(prepared_path) if prepared_path else OUT/'prepared'; meta = json.loads((prepared/'frame.json').read_text())
    op = Operator(prepared, meta); operator_checks = op.check_adjoint(); quotient = Quotient(op)
    Q = [positive(q) for q in op.Q]
    # Physical dual = y0 - y. Initialize from the previous quartic dual,
    # projected onto all new free-ideal equalities; this is only an iterate.
    y = quotient.project(quotient.y0+op.initial_dual())
    if restart:
        raw = np.load(restart)
        for i in range(len(Q)):
            previous = raw[f'Q_{i}']
            if previous.shape[0] > Q[i].shape[0] or previous.shape[0] != previous.shape[1]:
                raise ValueError('Restart does not embed in declared compact spans')
            Q[i] = np.zeros_like(Q[i]); Q[i][:previous.shape[0], :previous.shape[1]] = previous
        y = quotient.project(raw['y'])
    ATy = quotient.AT(y)
    Z = [positive(c-a) for c, a in zip(quotient.C, ATy)]
    def normal(v):
        return quotient.A(quotient.AT(v))
    operator = LinearOperator((len(quotient.rhs), len(quotient.rhs)), matvec=normal)
    preconditioner = None
    conditioning_receipt = None
    if normal_cache:
        from scipy.linalg import cho_solve
        from research.reconstruction_compression_20260914.inputs import sha
        cache = Path(normal_cache); conditioning_receipt = json.loads((cache/'receipt.json').read_text())
        if (sha(prepared/'frame.json') != conditioning_receipt['prepared_frame_sha256'] or
            sha(cache/'cholesky.npy') != conditioning_receipt['factor_sha256']):
            raise ValueError('Changed normal preconditioner')
        factor = np.load(cache/'cholesky.npy')
        inverse = lambda v: cho_solve((factor, True), v, check_finite=False)
        gv = inverse(quotient.v); denominator = float(quotient.v@gv)
        def conditioned(v):
            v = quotient.project(v)
            answer = inverse(v)-gv*(gv@v)/denominator
            return quotient.project(answer)
        preconditioner = LinearOperator(operator.shape, matvec=conditioned)
    # Diagonal scaling is applied in the outer coefficient maps. The exact
    # ideal projection does not commute with arbitrary diagonal left scaling,
    # so this CG uses no extra diagonal preconditioner.
    history = []; best = {'score': -1e100}; last_log = time.monotonic()
    def checkpoint(iteration, cg_status, cg_iterations):
        x = quotient.recover(Q)
        residual = (op.A(Q)+op.free@x-op.rhs)/op.scale
        r1 = float(np.sum(abs(residual)))
        dual = [c-a-z for c, a, z in zip(quotient.C, quotient.AT(y), Z)]
        dual_norm = float(np.sqrt(sum(np.sum(a*a) for a in dual)))
        rec = {'iteration': iteration, 'seconds': time.monotonic()-started, 'b': float(x[0]),
               'selected_row_l1': r1, 'projected_primal_l2': float(np.linalg.norm(quotient.A(Q)-quotient.rhs)),
               'dual_residual_l2': dual_norm, 'mu': mu, 'cg_status': cg_status, 'cg_iterations': cg_iterations,
               'untrusted_dual_ceiling': float(quotient.offset-quotient.rhs@y),
               'untrusted_objective_gap': float(quotient.offset-quotient.rhs@y-x[0])}
        score = x[0]-8*r1
        if score > best['score']:
            best.update(score=score, Q=[q.copy() for q in Q], x=x.copy(), record=rec)
        history.append(rec); print(json.dumps(rec), flush=True); dump(folder/'history.json', history)
        np.savez_compressed(folder/'checkpoint.npz', y=y, x=x, **{f'Q_{i}': q for i, q in enumerate(Q)})
        return rec
    details = {'kind': 'eliminated_full_number_quotient', 'gram_entries': sum(q.size for q in Q),
               'coefficient_rows_before_numerical_elimination': len(op.rhs), 'setup_seconds': time.monotonic()-started,
               'all_compact_span_cross_terms_retained': True, 'independent_partner_Gram_blocks': True,
               'many_body_states_enumerated': 0, 'full_Gram_teacher_used': False, 'additional_MPS_moments_used': False,
               'quotient_conditioning': quotient.record, 'matrix_free_adjoint_check': operator_checks,
               'fixed_mu': fixed_mu, 'prepared': str(prepared),
               'exact_number_dressed_linear_completion': meta.get('exact_number_dressed_linear_completion'),
               'normal_preconditioner': None if conditioning_receipt is None else {
                   k: v for k, v in conditioning_receipt.items() if k != 'blocks'}}
    dump(folder/'construction.json', details); print(json.dumps(details), flush=True)
    initial = checkpoint(0, 0, 0)
    best['score'] = -1e100  # Keep the inherited control separately; select an actual new iterate.
    solve_start = time.monotonic()
    for iteration in range(10000):
        if time.monotonic()-solve_start >= seconds:
            break
        residual = quotient.rhs-quotient.A(Q)
        rhs_y = quotient.A([c-z for c, z in zip(quotient.C, Z)])+residual/mu
        counter = [0]
        def count(_): counter[0] += 1
        y, status = cg(operator, rhs_y, x0=y, rtol=1e-7, atol=1e-12, maxiter=max_cg, M=preconditioner, callback=count)
        y = quotient.project(y)
        ATy = quotient.AT(y)
        Znew = []; Qnew = []
        for q, c, a in zip(Q, quotient.C, ATy):
            W = c-a-q/mu
            z = positive(W)
            qnew = q+mu*(a+z-c)
            Qnew.append((qnew+qnew.T)/2); Znew.append(z)
        Q, Z = Qnew, Znew
        primal_norm = float(np.linalg.norm(quotient.A(Q)-quotient.rhs))
        dual_norm = float(np.sqrt(sum(np.sum((c-a-z)**2) for c, a, z in zip(quotient.C, ATy, Z))))
        if iteration == 0 or time.monotonic()-last_log >= 10:
            rec = checkpoint(iteration+1, int(status), counter[0]); last_log = time.monotonic()
        if not fixed_mu and (iteration+1) % 25 == 0:
            if primal_norm > 5*dual_norm: mu = min(mu*1.5, 1000)
            elif dual_norm > 5*primal_norm: mu = max(mu/1.5, 1e-6)
        if primal_norm < 1e-9 and dual_norm < 1e-9:
            break
    checkpoint('final', int(status), counter[0])
    ex = export(op, best['Q'], best['x'], folder/'best')
    last = export(op, Q, quotient.recover(Q), folder/'last')
    details.update(best=ex, last=last, initial_control=initial, total_seconds=time.monotonic()-started,
                   status='requires_independent_exact_replay')
    dump(folder/'discovery.json', details)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('tag'); p.add_argument('--seconds', type=int, default=300)
    p.add_argument('--max-cg', type=int, default=100); p.add_argument('--mu', type=float, default=.1); p.add_argument('--restart')
    p.add_argument('--normal-cache')
    p.add_argument('--fixed-mu', action='store_true')
    p.add_argument('--prepared')
    a = p.parse_args(); run(a.tag, a.seconds, a.max_cg, a.mu, a.restart, a.normal_cache, a.fixed_mu, a.prepared)
