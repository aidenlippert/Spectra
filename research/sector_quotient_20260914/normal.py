"""Measured numerical conditioning: a bounded dense coefficient normal matrix.

The dense object is on polynomial coefficient equations, never N-particle
determinants. It is a costly reusable numerical preconditioner, not a compact
certificate. Partner cancellation halves the largest matrix multiplications.
"""
import json
import time
import numpy as np
from scipy import linalg
from research.sector_quotient_20260914.search import Operator
from research.sector_quotient_20260914.eliminated import Quotient
from research.sector_quotient_20260914.budget import OUT
from research.reconstruction_compression_20260914.inputs import dump, sha


def projected(M, V):
    active = np.flatnonzero(np.asarray(M.getnnz(axis=1)).ravel())
    M = M[active, :].tocsc(); n, r = V.shape; ii, jj = np.triu_indices(r)
    P = np.empty((len(active), len(ii)), dtype=np.float64)
    if n == r and np.array_equal(V, np.eye(n)):
        P[:] = (M[:, ii*n+jj]+M[:, jj*n+ii]).toarray()
        P[:, ii == jj] *= .5
    else:
        for start in range(0, len(ii), 24):
            i, j = ii[start:start+24], jj[start:start+24]
            X = np.einsum('ib,jb->ijb', V[:, i], V[:, j])
            X += np.einsum('ib,jb->ijb', V[:, j], V[:, i])*(i != j)[None, None, :]
            P[:, start:start+24] = M@X.reshape(n*n, -1)
    P[:, ii != jj] /= np.sqrt(2)
    return active, P


def add_rectangle(G, rows, cols, C, weight=1.):
    for start in range(0, len(rows), 192):
        rr = rows[start:start+192]; index = np.ix_(rr, cols)
        G[index] = G[index]+weight*C[start:start+192]


def run(tag='normal', prepared_path=None):
    if (OUT/'manifest.json').exists(): raise RuntimeError('Sealed campaign')
    started = time.monotonic(); folder = OUT/tag; folder.mkdir(exist_ok=False)
    from pathlib import Path
    prepared = Path(prepared_path) if prepared_path else OUT/'prepared'; meta = json.loads((prepared/'frame.json').read_text())
    op = Operator(prepared, meta); quotient = Quotient(op)
    if len(op.rhs) > 9000: raise ValueError('Dense numerical conditioning dimension budget')
    G = np.zeros((len(op.rhs), len(op.rhs)), order='F'); position = 0; blocks = []
    selected = np.load(prepared/'selected.npy')
    sextic = np.array([len(meta['rows'][i]) == 6 for i in selected])
    for k, pair in enumerate(meta['pairs']):
        t = time.monotonic(); V = op.V[position]; M = op.M[position]
        active, P = projected(M, V)
        if len(pair['members']) == 1:
            add_rectangle(G, active, active, P@P.T)
            position += 1
        else:
            n = V.shape[0]; order = np.array(pair['adjoint_order'])
            # Odd adjoint cancellation transposes the Gram pair indices.
            # On symmetric Gram inputs this leaves the operator unchanged.
            partner = op.M[position+1][:, (order[None, :]*n+order[:, None]).ravel()]
            paired = M+partner; paired.eliminate_zeros()
            if paired[sextic, :].nnz:
                raise AssertionError('Paired sextic map failed exact floating cancellation')
            low, L = projected(paired, V)
            add_rectangle(G, active, active, P@P.T, 2.)
            C = P@L.T
            add_rectangle(G, active, low, C, -1.)
            add_rectangle(G, low, active, C.T, -1.)
            add_rectangle(G, low, low, L@L.T)
            position += 2
        rec = {'pair': k, 'active_rows': len(active), 'symmetric_Gram_coordinates': P.shape[1],
               'seconds': time.monotonic()-t, 'elapsed_seconds': time.monotonic()-started}
        blocks.append(rec)
        if k >= 34 or k % 8 == 0: print(json.dumps(rec), flush=True)
        del P
    rng = np.random.default_rng(991)
    errors = []
    for _ in range(4):
        y = rng.normal(size=len(op.rhs)); direct = op.A(op.AT(y))
        errors.append(float(np.linalg.norm(G@y-direct)/np.linalg.norm(direct)))
    if max(errors) > 1e-9: raise AssertionError('Dense coefficient normal matrix does not match matrix-free operator')
    U = quotient.U
    G -= U@(U.T@G)
    G -= (G@U)@U.T
    G += U@U.T
    G = np.asarray((G+G.T)*.5, order='F')
    regularizer = 1e-10; G[np.diag_indices(len(G))] += regularizer
    print(json.dumps({'stage': 'factor', 'seconds': time.monotonic()-started, 'normal_probe_relative_errors': errors}), flush=True)
    C = linalg.cholesky(G, lower=True, overwrite_a=True, check_finite=False)
    np.save(folder/'cholesky.npy', C)
    np.save(folder/'ideal_basis.npy', U)
    record = {'dimension': len(op.rhs), 'dense_bytes': G.nbytes, 'diagonal_regularizer': regularizer,
              'normal_probe_relative_errors': errors, 'blocks': blocks, 'seconds': time.monotonic()-started,
              'seed_sha256': meta['seed_sha256'], 'prepared_frame_sha256': sha(prepared/'frame.json'),
              'factor_sha256': sha(folder/'cholesky.npy'), 'ideal_basis_sha256': sha(folder/'ideal_basis.npy'),
              'many_body_states_enumerated': 0, 'role': 'Numerical conditioning only, all final residuals require exact replay'}
    dump(folder/'receipt.json', record); print(json.dumps({k: v for k, v in record.items() if k != 'blocks'}), flush=True)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('tag', nargs='?', default='normal'); p.add_argument('--prepared')
    a = p.parse_args(); run(a.tag, a.prepared)
