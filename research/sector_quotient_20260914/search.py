"""Compact full-quotient molecular search without a dense projected CAR table."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsqr
from experiments.marginal_symbolic import decode, encode, add, scale, product
from research.collective_completion_20260914.spin_replay import setup, alpha_shift
from research.collective_completion_20260914.spin_screen import spin_squared, ladder_ideal
from research.reconstruction_compression_20260914.inputs import dump, sha
from research.sector_quotient_20260914.budget import ROOT, OUT
from research.sector_quotient_20260914.boundary import solve, positive


class Operator:
    def __init__(self, prepared, meta, paired=False):
        self.prepared = Path(prepared)
        self.meta = meta
        self.free = sparse.load_npz(self.prepared/'free.npz')
        self.rhs = np.load(self.prepared/'rhs.npy')
        seed = np.load(self.prepared/'seed.npz')
        self.V = []; self.M = []; self.Q = []; self.members = []; self.ids = []
        norms = np.asarray(self.free.power(2).sum(axis=1)).ravel()
        for k, pair in enumerate(meta['pairs']):
            V = seed[f'V_{k}']
            matrices = []
            for j in pair['members']:
                M = sparse.load_npz(self.prepared/f'map_{j}.npz')
                matrices.append(M)
            if paired and len(matrices) == 2:
                order = np.array(pair['adjoint_order']); size = V.shape[0]
                M = matrices[0]+matrices[1][:, (order[:, None]*size+order[None, :]).ravel()]
                M.eliminate_zeros()
                self.V.append(V); self.M.append(M); self.Q.append(seed[f'Q_{k}'])
                self.members.append(pair['members']); self.ids.append(k)
                norms += np.asarray(M.power(2).sum(axis=1)).ravel()
            else:
                for a, (j, M) in enumerate(zip(pair['members'], matrices)):
                    W = V
                    if a:
                        W = np.empty_like(V); W[pair['adjoint_order']] = V
                    self.V.append(W); self.M.append(M); self.Q.append(seed[f'Q_{k}'])
                    self.members.append([j]); self.ids.append(k)
                    norms += np.asarray(M.power(2).sum(axis=1)).ravel()
        self.scale = 1/np.sqrt(np.maximum(1, norms))
        self.M = [(sparse.diags(self.scale)@M).tocsc() for M in self.M]
        self.free = (sparse.diags(self.scale)@self.free).tocsc()
        self.rhs = self.scale*self.rhs
        old = seed['old_x']; nnumber = len(meta['number_basis']); nspin = len(meta['spin_basis'])
        self.x = np.zeros(self.free.shape[1])
        self.x[0] = old[0]; self.x[1:1+nspin] = old[1:1+nspin]
        self.x[1+nnumber:] = old[1+nspin:]
        # A weighted Frobenius diagonal upper bound is cheap and positive.
        # A deterministic matrix-free Hutchinson estimate below refines it.
        self.diagonal = np.asarray(self.free.power(2).sum(axis=1)).ravel()
        self.identity = [np.array_equal(V, np.eye(V.shape[0])) if V.shape[0] == V.shape[1] else False for V in self.V]

    def A(self, matrices):
        out = np.zeros(len(self.rhs))
        for M, V, Q, identity in zip(self.M, self.V, matrices, self.identity):
            physical = Q if identity else V@Q@V.T
            out += M@physical.ravel()
        return out

    def AT(self, y):
        out = []
        for M, V, identity in zip(self.M, self.V, self.identity):
            C = np.asarray(M.T@y).reshape(V.shape[0], V.shape[0])
            if not identity:
                C = V.T@C@V
            out.append((C+C.T)/2)
        return out

    def initial_dual(self):
        raw = np.load(self.meta['seed'])
        yold = raw['dual']
        T = sparse.load_npz(self.prepared/'twirl.npz')
        selected = np.load(self.prepared/'selected.npy')
        small = sum(len(w) <= 4 for w in self.meta['rows'])
        selected_small = selected[selected < small]
        inverse = lsqr(T[selected_small, :small].T, yold, atol=1e-12, btol=1e-12)
        y = np.zeros(len(selected)); y[:len(selected_small)] = -inverse[0]
        return y/self.scale

    def check_adjoint(self):
        rng = np.random.default_rng(1917)
        Q = [(lambda a: (a+a.T)/2)(rng.normal(size=q.shape)) for q in self.Q]
        y = rng.normal(size=len(self.rhs))
        left = self.A(Q)@y
        right = sum(np.sum(a*q) for a, q in zip(self.AT(y), Q))
        relative = abs(left-right)/max(1., abs(left), abs(right))
        if relative > 1e-10:
            raise AssertionError('Matrix-free Frobenius adjoint failed')
        # The estimator is a numerical preconditioner, not a physical moment.
        diagonal = np.zeros(len(y))
        for _ in range(12):
            noise = [(lambda a: (a+a.T)/np.sqrt(2))(rng.normal(size=q.shape)) for q in self.Q]
            # Off-diagonal entries have variance 1; diagonal entries variance 2.
            # Divide the latter by sqrt(2), matching a symmetric Frobenius basis.
            for q in noise:
                q[np.diag_indices(len(q))] /= np.sqrt(2)
            v = self.A(noise)
            diagonal += v*v/12
        self.diagonal += diagonal
        return {'relative_adjoint_error': relative, 'preconditioner_probes': 12}


def export(operator, Qs, x, folder):
    meta = operator.meta
    data = json.loads(Path(meta['fixture']).read_text()); m, n = data['modes'], data['particles']
    h, hs, delta = setup(data)
    number_basis = [decode(p, m, 4) for p in meta['number_basis']]
    spin_basis = [decode(p, m, 2) for p in meta['spin_basis']]
    ladder_basis = [decode(p, m, 2) for p in meta['ladder_basis']]
    def poly(basis, values):
        return add(*(scale(p, F(round(float(value)*10**12), 10**12)) for p, value in zip(basis, values)))
    nn, ns = len(number_basis), len(spin_basis)
    X = poly(number_basis, x[1:1+nn]); Y = poly(spin_basis, x[1+nn:1+nn+ns])
    a = F(round(float(x[1+nn+ns])*10**12), 10**12)
    W = poly(ladder_basis, x[2+nn+ns:])
    blocks = []; min_eigenvalue = 0.
    for V, Q, members, kid in zip(operator.V, Qs, operator.members, operator.ids):
        ev, U = np.linalg.eigh((Q+Q.T)/2)
        min_eigenvalue = min(min_eigenvalue, float(ev[0]))
        R = np.sqrt(np.maximum(ev, 0))[:, None]*U.T
        factors = np.rint((R@V.T)*1e9).astype(np.int64)
        for a_index, i in enumerate(members):
            Z = factors
            if a_index:
                Z = np.empty_like(factors); Z[:, meta['pairs'][kid]['adjoint_order']] = factors
            Z = Z[np.any(Z, axis=1)]
            if len(Z):
                blocks.append({'name': meta['groups'][i]['name'], 'words': meta['groups'][i]['words'], 'factor': Z.tolist()})
    core = {'modes': m, 'particles': n, 'operator_degree': 3,
            'hamiltonian': encode(add(hs, scale(product(alpha_shift(m, n), Y), -1), scale(spin_squared(m), -a), scale(ladder_ideal(m, W), -1))),
            'b': str(F(round(float(x[0])*1e12), 10**12)), 'number_multiplier': encode(X), 'denominator': 10**9, 'blocks': blocks}
    cert = {'kind': 'spin_sector_sos_v1', 'modes': m, 'particles': n, 'hamiltonian': encode(h),
            'alpha_multiplier': encode(Y), 'core': core, 'magnetization': 0, 'singlet': True,
            'casimir_multiplier': str(a), 'spin_ladder_multiplier': encode(W), 'spin_twirl': True}
    folder.mkdir(exist_ok=False)
    with (folder/'certificate.json').open('x') as f:
        json.dump(cert, f, separators=(',', ':'))
    np.savez_compressed(folder/'raw.npz', x=x, **{f'Q_{k}': Q for k, Q in enumerate(Qs)})
    return {'proposed_b': float(x[0]), 'number_multiplier_quartic_terms': sum(len(t['word']) == 4 for t in core['number_multiplier']),
            'minimum_eigenvalue_before_export': min_eigenvalue, 'factor_rows': sum(len(b['factor']) for b in blocks),
            'certificate_bytes': (folder/'certificate.json').stat().st_size}


def run(tag, seconds=240, max_cg=60, mu=1., paired=False, restart=None):
    if (OUT/'manifest.json').exists():
        raise RuntimeError('Sealed campaign')
    total = time.monotonic(); prepared = OUT/'prepared'
    meta = json.loads((prepared/'frame.json').read_text())
    if sha(meta['fixture']) != meta['fixture_sha256'] or sha(meta['seed']) != meta['seed_sha256']:
        raise ValueError('Changed input')
    folder = OUT/'candidates'/tag; folder.mkdir(parents=True, exist_ok=False)
    op = Operator(prepared, meta, paired=paired)
    checks = op.check_adjoint()
    raw_residual = (op.A(op.Q)+op.free@op.x-op.rhs)/op.scale
    if np.max(abs(raw_residual)) > 1e-5:
        raise AssertionError('Paired seed lost in the enlarged formulation')
    details = {'kind': 'matrix_free_boundary_proposal', 'paired_control': paired,
               'coefficient_equations': len(op.rhs), 'free_variables': op.free.shape[1],
               'gram_dimensions': [q.shape[0] for q in op.Q], 'gram_entries': sum(q.size for q in op.Q),
               'sparse_physical_map_nonzeros': sum(M.nnz for M in op.M),
               'projected_dense_map_materialized': False, 'paired_seed_max_coefficient_error': float(max(abs(raw_residual))),
               'all_compact_span_cross_terms_retained': True, 'independent_partner_Gram_blocks': not paired,
               'quartic_number_basis_dimension': meta['quartic_number_basis_dimension'],
               'additional_MPS_moments_used': False, 'full_Gram_teacher_used': False, 'many_body_states_enumerated': 0,
               'setup_seconds': time.monotonic()-total, **checks}
    dump(folder/'construction.json', details); print(json.dumps(details), flush=True)
    best = {'score': -1e100, 'Q': None, 'x': None, 'record': None}
    history = []
    def callback(record, Q, x, y, Z):
        defect = (op.A(Q)+op.free@x-op.rhs)/op.scale
        # This is only candidate selection. Exact CAR replay alone accepts it.
        l1 = float(np.sum(abs(defect)))
        score = float(x[0])-8*l1
        record = {**record, 'selected_row_l1': l1, 'untrusted_selection_score': score}
        history.append(record); print(json.dumps(record), flush=True)
        if score > best['score']:
            best.update(score=score, Q=[q.copy() for q in Q], x=x.copy(), record=record)
        np.savez_compressed(folder/'checkpoint.npz', x=x, y=y, **{f'Q_{i}': q for i, q in enumerate(Q)}, **{f'Z_{i}': z for i, z in enumerate(Z)})
        dump(folder/'history.json', history)
    y = op.initial_dual()
    if restart:
        raw = np.load(restart); op.x = raw['x']; op.Q = [raw[f'Q_{i}'] for i in range(len(op.Q))]; y = raw['y']
    Q, x, y, Z, records = solve(op.A, op.AT, op.free, op.rhs, op.Q, op.x, seconds,
                               callback=callback, mu=mu, max_cg=max_cg, diagonal=op.diagonal, initial_y=y)
    callback({'iteration': 'final', 'seconds': time.monotonic()-total, 'b': float(x[0])}, Q, x, y, Z)
    best_export = export(op, best['Q'], best['x'], folder/'best')
    last_export = export(op, Q, x, folder/'last')
    details.update(best=best_export, last=last_export, total_seconds=time.monotonic()-total,
                   status='requires_independent_exact_replay', maximum_cg_iterations=max_cg, requested_solve_seconds=seconds)
    dump(folder/'discovery.json', details)
    print(json.dumps({k: v for k, v in details.items() if k != 'gram_dimensions'}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('tag'); p.add_argument('--seconds', type=int, default=240)
    p.add_argument('--max-cg', type=int, default=60); p.add_argument('--mu', type=float, default=1.)
    p.add_argument('--paired', action='store_true'); p.add_argument('--restart')
    a = p.parse_args(); run(a.tag, a.seconds, a.max_cg, a.mu, a.paired, a.restart)
