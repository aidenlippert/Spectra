"""Optimize the actual lifted coefficient-L1 bound on a declared SOS subcone.

Ideals-only mode keeps every Gram matrix fixed. Block mode permits nonnegative
rescaling of those same PSD blocks jointly with all independent sector ideals.
Neither mode asserts an optimum of the full Gram family. Exact replay accepts.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu
from research.acceptance_channels_20260915.campaign import ROOT
from research.gpu_acceleration_20260915.solve import Operator
from research.sector_quotient_20260914.search import export


def residual_lift(prepared, scale):
    T = sparse.load_npz(prepared / 'twirl.npz').tocsc()
    selected = np.load(prepared / 'selected.npy')
    J = T[selected][:, selected].tocsr()
    count, labels = connected_components(J + J.T, directed=False)
    order = np.argsort(labels, kind='stable')
    sizes = np.bincount(labels)
    boundaries = np.r_[0, sizes.cumsum()]
    if sizes.max(initial=0) > 64:
        raise ValueError('Residual lift exceeds the declared small-block envelope')
    rows, cols, values = [], [], []
    cache = {}
    for begin, end in zip(boundaries[:-1], boundaries[1:]):
        ids = order[begin:end]
        block = J[ids][:, ids].toarray()
        key = (block.shape, block.tobytes())
        if key not in cache:
            cache[key] = np.linalg.inv(block)
        inverse = cache[key]
        i, j = np.nonzero(inverse)
        rows.extend(ids[i]); cols.extend(ids[j]); values.extend(inverse[i, j])
    inverse = sparse.csc_matrix((values, (rows, cols)), shape=J.shape)
    defect = J @ inverse - sparse.eye(J.shape[0], format='csc')
    if max(abs(defect.data), default=0) > 1e-12:
        raise ValueError('Numerical residual lift failed its inverse check')
    K = (T[:, selected] @ inverse @ sparse.diags(1 / scale)).tocsr()
    rng = np.random.default_rng(714)
    probe = rng.normal(size=len(scale))
    expected = T[:, selected] @ splu(J.tocsc()).solve(probe / scale)
    if np.max(abs(K @ probe - expected)) > 1e-9:
        raise ValueError('Lift differs from the existing accepting-score reconstruction')
    return K, {
        'spin_components': int(count), 'largest_spin_component': int(sizes.max()),
        'distinct_numerical_inverse_patterns': len(cache),
        'lift_shape': list(K.shape), 'lift_nonzeros': K.nnz,
        'dense_global_inverse_built': False,
    }


def l1_fit(L, residual, weights, gain, lower_bounds, seconds):
    """Minimize ||residual-L v||_(weights,1)-gain.v; numerical proposal only."""
    L = sparse.csr_matrix(L)
    residual, weights, gain = map(np.asarray, (residual, weights, gain))
    if L.shape != (len(residual), len(gain)) or weights.shape != residual.shape:
        raise ValueError('Coefficient-fit dimensions disagree')
    if not all(np.isfinite(v).all() for v in (L.data, residual, weights, gain)) or np.any(weights <= 0):
        raise ValueError('Finite coefficients and strictly positive weights required')
    active = np.flatnonzero(np.diff(L.indptr))
    constant = float(weights @ abs(residual) - weights[active] @ abs(residual[active]))
    if len(active) == 0:
        raise ValueError('No variable affects the residual')
    # Work in mHa-sized increments to keep LP feasibility tolerances meaningful.
    unit = .001
    matrix = L[active]
    rhs = residual[active] / unit
    I = sparse.eye(len(active), format='csr')
    inequalities = sparse.vstack((sparse.hstack((-matrix, -I)), sparse.hstack((matrix, -I))), format='csc')
    objective = np.r_[-gain, weights[active]]
    bounds = [(None if v is None else v / unit, None) for v in lower_bounds]
    bounds += [(0, None)] * len(active)
    result = linprog(objective, A_ub=inequalities, b_ub=np.r_[-rhs, rhs],
        bounds=bounds, method='highs-ds', options={'time_limit': seconds,
        'dual_feasibility_tolerance': 1e-9, 'primal_feasibility_tolerance': 1e-9})
    details = {'lp_status': int(result.status), 'lp_message': result.message,
        'active_rows': len(active), 'fixed_residual_allowance_Ha': constant,
        'variables_without_epigraph': L.shape[1], 'inequality_nonzeros': inequalities.nnz,
        'iterations': int(result.nit), 'exact_family_obstruction': False}
    if result.x is None:
        return None, details
    step = unit * result.x[:L.shape[1]]
    if any(bound is not None and value < bound-1e-8 for value, bound in zip(step, lower_bounds)):
        raise ValueError('Numerical coefficient fit violated a positivity bound')
    details['objective_Ha'] = float(weights @ abs(residual-L @ step) - gain @ step)
    return step, details


def run(case, source, tag, mode, seconds):
    start = time.monotonic()
    out = case / tag
    out.mkdir(exist_ok=False)
    sys.path.insert(0, str(ROOT / '.venv-interacting-libs'))
    from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient
    op = Operator(case)
    fixture = case / 'fixture.json'
    if hashlib.sha256(fixture.read_bytes()).hexdigest() != op.meta['fixture_sha256']:
        raise ValueError('Hamiltonian changed')
    op.meta['fixture'] = str(fixture)
    raw = np.load(source)
    Q = [raw[f'Q_{i}'] for i in range(len(op.Q))]
    x = raw['x'].copy()
    if any(q.shape != expected.shape for q, expected in zip(Q, op.Q)):
        raise ValueError('Source Gram dimensions changed')
    quo = SparseQuotient(op)
    K, lift_details = residual_lift(op.prepared, op.scale)
    weights = np.load(op.prepared / 'weights.npy')
    residual = K @ (op.rhs - op.A(Q) - op.free @ x)
    ids = np.r_[0, quo.pivots]
    columns = [K @ op.free[:, ids]]
    gram_count = len(Q) if mode == 'blocks' else 0
    if gram_count:
        # Each column is one whole PSD block; retain its internal cross terms.
        columns += [sparse.hstack([sparse.csc_matrix((K @ (M @ q.ravel()))[:, None])
            for M, q in zip(op.M, Q)], format='csc')]
    L = sparse.hstack(columns, format='csr')
    gain = np.zeros(L.shape[1]); gain[0] = 1
    lower = [None] * len(ids) + [-1.] * gram_count
    baseline_lower = float(x[0] - weights @ abs(residual))
    setup_seconds = time.monotonic()-start
    print(json.dumps({'stage': 'fit', 'mode': mode, 'setup_seconds': setup_seconds,
        'baseline_predicted_lower_Ha': baseline_lower, **lift_details}), flush=True)
    step, details = l1_fit(L, residual, weights, gain, lower, seconds)
    if step is not None:
        trial_x = x.copy(); trial_x[ids] += step[:len(ids)]
        factors = np.maximum(0., 1+step[len(ids):]) if gram_count else np.ones(len(Q))
        trial_Q = [q * factor for q, factor in zip(Q, factors)]
        new_residual = K @ (op.rhs-op.A(trial_Q)-op.free @ trial_x)
        candidate_lower = float(trial_x[0]-weights @ abs(new_residual))
        if candidate_lower > baseline_lower:
            Q, x, residual = trial_Q, trial_x, new_residual
    eta = float(weights @ abs(residual))
    U = float(Fraction(op.meta['upper_Ha']))
    delta = float(Fraction(op.meta['original_H_spin_defect_Ha']))
    export_receipt = export(op, Q, x, out / 'export')
    details.update(kind='accepted_l1_fit_on_fixed_gram_subcone', mode=mode,
        source=str(source), setup_seconds=setup_seconds, total_seconds=time.monotonic()-start,
        predicted_residual_Ha=eta, predicted_width_mHa=1000*(U-x[0]+eta+delta),
        predicted_improvement_mHa=1000*(x[0]-eta-baseline_lower),
        full_Gram_family_optimized=False, exact_acceptance_required=True,
        lift=lift_details, export=export_receipt)
    (out / 'discovery.json').write_text(json.dumps(details, indent=2)+'\n')
    print(json.dumps(details), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path); parser.add_argument('source', type=Path)
    parser.add_argument('tag'); parser.add_argument('--mode', choices=('ideals', 'blocks'), default='ideals')
    parser.add_argument('--seconds', type=float, default=180)
    args = parser.parse_args()
    run(args.case.resolve(), args.source.resolve(), args.tag, args.mode, args.seconds)
