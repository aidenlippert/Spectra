"""Same-family Douglas--Rachford feasibility at a declared scalar target.

This is a numerical construction experiment, not a family-limit oracle.
Every export uses the original exact checker and all its residual allowances.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator, cg, splu
from research.gpu_acceleration_20260915.solve import Operator
from research.gpu_acceleration_20260915.kernel_bench import cpu_project
from research.sector_quotient_20260914.search import export


def step(A, AT, F, target, W, z, solve, project):
    Q = project(W)
    reflected = [2*q-w for q, w in zip(Q, W)]
    direction = solve(target-A(reflected)-F@z)
    return [q+a for q, a in zip(Q, AT(direction))], z+F.T@direction


def run(case, source, tag, seconds, gap_mha):
    start = time.monotonic()
    out = case / tag
    out.mkdir(exist_ok=False)
    op = Operator(case)
    fixture = case / 'fixture.json'
    if hashlib.sha256(fixture.read_bytes()).hexdigest() != op.meta['fixture_sha256']:
        raise ValueError('Hamiltonian changed')
    op.meta['fixture'] = str(fixture)
    raw = np.load(source)
    W = [raw[f'Q_{i}'].copy() for i in range(len(op.Q))]
    if any(w.shape != q.shape for w, q in zip(W, op.Q)):
        raise ValueError('Source Gram dimensions changed')
    U = float(Fraction(op.meta['upper_Ha']))
    delta = float(Fraction(op.meta['original_H_spin_defect_Ha']))
    b = U-gap_mha/1000
    F = op.free[:, 1:].tocsr()
    f0 = op.free[:, 0]
    target = op.rhs-b*f0.toarray().ravel()
    G = (op.G-f0@f0.T).tocsc()
    pre = LinearOperator(G.shape, matvec=splu(G+sparse.eye(G.shape[0], format='csc')*1e-11).solve)
    last_cg = [0]
    def solve(rhs):
        value, status = cg(G, rhs, M=pre, rtol=1e-11, atol=1e-14, maxiter=80)
        if status < 0 or not np.isfinite(value).all():
            raise ValueError('Affine projection failed')
        last_cg[0] = int(status)
        return value
    rng = np.random.default_rng(19)
    probe = rng.normal(size=len(target))
    error = np.linalg.norm(G@probe-op.A(op.AT(probe))-F@(F.T@probe))/np.linalg.norm(G@probe)
    if error > 1e-10:
        raise ValueError('Fixed-scalar normal map disagrees with coefficient maps')
    T = sparse.load_npz(op.prepared/'twirl.npz')
    selected = np.load(op.prepared/'selected.npy')
    Tc = T[:, selected].tocsc()
    inverse = splu(T[selected][:, selected].tocsc())
    weights = np.load(op.prepared/'weights.npy')
    project = lambda blocks: cpu_project(blocks, 'evd')
    z = raw['x'][1:].copy()
    best = {'score': -float('inf')}
    records = []
    def record(iteration):
        Q = project(W)
        x = np.r_[b, z]
        r = (op.rhs-op.A(Q)-op.free@x)/op.scale
        eta = float(weights@abs(Tc@inverse.solve(r)))
        score = b-eta-delta
        value = {'iteration': iteration, 'seconds': time.monotonic()-start,
            'b': b, 'scalar_gap_mHa': gap_mha, 'predicted_residual_Ha': eta,
            'unverified_width_mHa': 1000*(U-score), 'cg_status': last_cg[0]}
        records.append(value)
        if score > best['score']:
            best.update(score=score, Q=[q.copy() for q in Q], x=x.copy(), record=value)
        (out/'history.json').write_text(json.dumps(records, indent=2)+'\n')
        print(json.dumps(value), flush=True)
        return value
    record(0)
    setup_seconds = time.monotonic()-start
    iteration = 0
    while time.monotonic()-start < seconds:
        W, z = step(op.A, op.AT, F, target, W, z, solve, project)
        iteration += 1
        if iteration == 1 or iteration % 100 == 0:
            value = record(iteration)
            if value['unverified_width_mHa'] < 1.5:
                break
    record('final')
    result = export(op, best['Q'], best['x'], out/'export')
    np.savez_compressed(out/'checkpoint.npz', z=z, **{f'W_{i}': w for i, w in enumerate(W)})
    details = {'kind': 'fixed_scalar_affine_PSD_feasibility', 'source': str(source),
        'best': best['record'], 'setup_seconds': setup_seconds,
        'total_seconds': time.monotonic()-start, 'normal_relative_error': float(error),
        'same_Gram_family': True, 'fixed_scalar_gap_mHa': gap_mha,
        'family_obstruction_proved': False, 'exact_acceptance_required': True, 'export': result}
    (out/'discovery.json').write_text(json.dumps(details, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path); parser.add_argument('source', type=Path)
    parser.add_argument('tag'); parser.add_argument('--seconds', type=float, default=300)
    parser.add_argument('--gap-mha', type=float, default=1.)
    args = parser.parse_args()
    run(args.case.resolve(), args.source.resolve(), args.tag, args.seconds, args.gap_mha)
