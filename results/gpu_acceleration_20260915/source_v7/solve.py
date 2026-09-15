"""Isolated acceleration experiment; acceptance remains the original exact checker.

Derived from the frozen transfer solver. Backend-only runs preserve its equations
and tolerances. Fixed iteration schedules make CPU/GPU timing comparable.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import LinearOperator, cg, splu
from research.sector_quotient_20260914.eliminated import Quotient
from research.sector_quotient_20260914.search import export
from research.transfer_solver_20260915.budget import dump


class Operator:
    def __init__(self, case):
        self.prepared = case/'prepared'
        self.meta = json.loads((self.prepared/'frame.json').read_text())
        z = np.load(self.prepared/'bases.npz')
        self.V = [z[f'V_{i}'] for i in range(len(self.meta['blocks']))]
        self.members = [[b['physical_group']] for b in self.meta['blocks']]
        self.ids = [0]*len(self.V)
        self.Q = [np.zeros((V.shape[1], V.shape[1])) for V in self.V]
        self.M = [sparse.load_npz(self.prepared/f'map_{i}.npz').tocsr() for i in range(len(self.V))]
        self.MT = [M.T.tocsr() for M in self.M]
        self.free = sparse.load_npz(self.prepared/'free.npz').tocsc()
        self.rhs = np.load(self.prepared/'rhs.npy')
        self.scale = np.load(self.prepared/'scale.npy')
        self.x = np.zeros(self.free.shape[1])
        self.G = sparse.load_npz(self.prepared/'normal.npz').tocsc()

    def A(self, Q):
        return sum((M@q.ravel() for M, q in zip(self.M, Q)), np.zeros(len(self.rhs)))

    def AT(self, y):
        return [np.asarray(M@y).reshape(V.shape[1], V.shape[1]) for M, V in zip(self.MT, self.V)]


def run(case, tag, seconds, mu, restart, backend='cpu_evr', iterations=None,
        record_every=100, switch_iteration=None, adaptive=False):
    started = time.monotonic()
    out = case/tag
    out.mkdir(exist_ok=False)
    op = Operator(case)
    # Only a relocated, copied fixture path changes; the hash must still match.
    import hashlib
    fixture = case/'fixture.json'
    if hashlib.sha256(fixture.read_bytes()).hexdigest() != op.meta['fixture_sha256']:
        raise ValueError('Changed Hamiltonian input')
    op.meta['fixture'] = str(fixture)
    Q, x = op.Q, op.x
    y = -np.load(op.prepared/'physical_dual.npy')
    if restart:
        z = np.load(restart)
        Q = [z[f'Q_{i}'] for i in range(len(Q))]
        x, y = z['x'], z['y']
        if any(q.shape != old.shape for q, old in zip(Q, op.Q)):
            raise ValueError('Restart representation mismatch')
    c = np.zeros(op.free.shape[1])
    c[0] = -1.
    from research.gpu_acceleration_20260915.kernel_bench import cpu_project, GPUProject, HybridProject
    if backend == 'hybrid':
        projector = HybridProject([q.shape for q in Q])
    elif backend.startswith('gpu'):
        projector = GPUProject([q.shape for q in Q], batched=backend != 'gpu_individual', mixed=backend == 'gpu_mixed')
    else:
        projector = lambda blocks: cpu_project(blocks, backend.removeprefix('cpu_'))
    Z = projector([-v for v in op.AT(y)])
    lu = splu(op.G+sparse.eye(len(op.rhs), format='csc')*1e-11)
    pre = LinearOperator(op.G.shape, matvec=lu.solve)
    rng = np.random.default_rng(51)
    probe = rng.normal(size=len(op.rhs))
    normal_error = np.linalg.norm(op.G@probe-op.A(op.AT(probe))-op.free@(op.free.T@probe))/np.linalg.norm(op.G@probe)
    if normal_error > 1e-10:
        raise ValueError('Normal system does not match declared coefficient maps')
    quo = Quotient(op)
    T = sparse.load_npz(op.prepared/'twirl.npz')
    selected = np.load(op.prepared/'selected.npy')
    Tc = T[:, selected].tocsc()
    residual_inverse = splu(T[selected][:, selected].tocsc())
    weights = np.load(op.prepared/'weights.npy')
    U = float(F(op.meta['upper_Ha']))
    delta = float(F(op.meta['original_H_spin_defect_Ha']))
    best = {'score': -float('inf')}
    history = []
    setup_seconds = time.monotonic()-started
    kernel_seconds = 0.
    def record(iteration, cg_status=0):
        xx = quo.recover(Q)
        r = (op.A(Q)+op.free@xx-op.rhs)/op.scale
        full = -Tc@residual_inverse.solve(r)
        if np.max(abs(full[selected]+r)) > 1e-7:
            raise ValueError('Residual reconstruction lost coefficient constraints')
        eta = float(weights@abs(full))
        score = float(xx[0])-eta-delta
        rec = {'iteration': iteration, 'seconds': time.monotonic()-started, 'mu': mu,
               'b': float(xx[0]), 'raw_width_mHa': 1000*(U-xx[0]),
               'predicted_coefficient_l1_Ha': eta, 'unverified_width_mHa': 1000*(U-score),
               'primal_l2': float(np.linalg.norm(op.A(Q)+op.free@x-op.rhs)),
               'dual_l2': float(np.sqrt(sum(np.sum((a+z)**2) for a, z in zip(op.AT(y), Z))
                                          + np.linalg.norm(op.free.T@y-c)**2)),
               'minimum_dual_eigenvalue': min(float(linalg.eigvalsh(-a, subset_by_index=(0, 0), check_finite=False)[0]) for a in op.AT(y)),
               'cg_status': int(cg_status), 'kernel_seconds': kernel_seconds}
        history.append(rec)
        (out/'history.json').write_text(json.dumps(history, indent=2)+'\n')
        np.savez_compressed(out/'checkpoint.npz', x=x, y=y, **{f'Q_{i}': q for i, q in enumerate(Q)})
        if score > best['score']:
            best.update(score=score, record=rec, Q=[q.copy() for q in Q], x=xx.copy())
        print(json.dumps(rec), flush=True)
        return rec
    record(0)
    last = time.monotonic()
    cg_status = 0
    switched = False
    completed = 0
    for it in range(1, 100000):
        if time.monotonic()-started >= seconds or (iterations is not None and it > iterations):
            break
        if not switched and ((switch_iteration is not None and it == switch_iteration) or
                (adaptive and history[-1]['unverified_width_mHa'] < 10 and
                 history[-1]['primal_l2'] < 1e-4)):
            mu = .03
            Z = projector([-v for v in op.AT(y)])
            switched = True
        step_started = time.monotonic()
        r = op.rhs-op.A(Q)-op.free@x
        target = op.free@c-op.A(Z)+r/mu
        y, cg_status = cg(op.G, target, x0=y, M=pre, rtol=1e-10, atol=1e-13, maxiter=40)
        if not np.isfinite(y).all():
            raise ValueError('Nonfinite boundary-point iterate')
        W = [q+mu*a for q, a in zip(Q, op.AT(y))]
        Qn = projector(W)
        Zn = [(qp-w)/mu for qp, w in zip(Qn, W)]
        x = x+mu*(op.free.T@y-c)
        Q, Z = Qn, Zn
        completed = it
        kernel_seconds += time.monotonic()-step_started
        if it == 1 or (it % record_every == 0 if iterations is not None or adaptive else time.monotonic()-last >= 10):
            rec = record(it, cg_status)
            last = time.monotonic()
            if iterations is None and rec['unverified_width_mHa'] < 1.3 and rec['primal_l2'] < 5e-8:
                break
    record('final', cg_status)
    result = export(op, best['Q'], best['x'], out/'export')
    dump(out/'discovery.json', {'kind': 'input_derived_sparse_boundary_point_v1', 'best': best['record'],
                              'export': result, 'seconds': time.monotonic()-started,
                              'normal_map_relative_error': float(normal_error),
                              'ideal_projection': quo.record, 'restart': str(restart) if restart else None,
                              'initial_old_checkpoint_used': False,
                              'Gram_entries': sum(q.size for q in Q),
                              'backend': backend, 'setup_seconds': setup_seconds,
                              'kernel_seconds': kernel_seconds, 'iterations_completed': completed,
                              'fixed_iterations_requested': iterations,
                              'switch_iteration': switch_iteration,
                              'adaptive_refinement': adaptive,
                              'status': 'requires_original_exact_checker'})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('tag')
    p.add_argument('--seconds', type=float, required=True)
    p.add_argument('--mu', type=float, required=True)
    p.add_argument('--restart', type=Path)
    p.add_argument('--backend', choices=('cpu_evr', 'cpu_evd', 'cpu_positive', 'gpu_batched', 'gpu_individual', 'gpu_mixed', 'hybrid'), default='cpu_evr')
    p.add_argument('--iterations', type=int)
    p.add_argument('--record-every', type=int, default=100)
    p.add_argument('--switch-iteration', type=int)
    p.add_argument('--adaptive', action='store_true')
    a = p.parse_args()
    run(a.case.resolve(), a.tag, a.seconds, a.mu, a.restart, a.backend,
        a.iterations, a.record_every, a.switch_iteration, a.adaptive)
