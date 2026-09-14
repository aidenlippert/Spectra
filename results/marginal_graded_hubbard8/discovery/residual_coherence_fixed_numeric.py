"""Bounded two-coordinate spectral probe; no production energy acceptance."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib
import json
import math
import sys
import time
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

from residual_coherence_obstruction import ROOT, BASE, OLD_MODULES, pure_coherence
from spin_word_numeric import prepare, legacy_terms
from experiments.marginal_spin_word_telescope import local_value
from experiments.marginal_hopping_telescope import projected_matrix


def action_from_matrix(matrix, denominator):
    result = [{} for _ in range(4096)]
    for (r, c), value in matrix.items():
        result[c][r] = F(value, denominator)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    certificate = json.loads(args.seed.read_text())
    candidates_path = BASE.parent / 'joint_projector/signed_density/symmetry_diagonal_candidates.json'
    candidates = json.loads(candidates_path.read_text())['candidates']
    shapes = [{int(s): v for s, v in candidates[i]['diagonal'].items()} for i in (1, 2, 4, 5, 6, 7, 31, 11, 20)]
    x, physical, blocks, *_ = prepare(certificate, shapes)
    old_actions, coefficients = [], []
    for name in OLD_MODULES:
        module = importlib.import_module('experiments.marginal_' + name)
        field = 'coherent_projector' if name == 'coherent_projector_telescope' else name
        for label in getattr(module, 'LABELS', ('single',)):
            old_actions.append(module.actions({label: 1}) if hasattr(module, 'LABELS') else module.actions())
            coefficients.append(F(certificate[field].get(label, 0)) if hasattr(module, 'LABELS') else F(certificate[field]))
    terms = legacy_terms(certificate)
    for label, value in certificate['spin_word_telescope'].items():
        terms[label] += F(value)
    alpha, beta = F(certificate['penalty']), F(certificate['joint']['penalty'])
    theta_h = F(certificate['projector_sum_ceiling']) / certificate['windows']
    theta_j = F(certificate['joint']['projector_sum_ceiling']) / certificate['joint']['windows']
    penalty_cost = alpha * theta_h + beta * theta_j
    new_actions, receipts = [], []
    for case in ('W_zero', 'W_plus_1'):
        path = BASE / case / 'telescope_replay.json'
        receipt = json.loads(path.read_text())
        if receipt.get('accepted') is not True:
            raise ValueError('Accepted new telescope required')
        for name, digest in receipt['source_sha256'].items():
            if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
                raise ValueError('Stale telescope source')
        vector = {int(s): a for s, a in receipt['five_site_vector'].items()}
        denominator, _, matrix, _ = pure_coherence(vector)
        new_actions.append(action_from_matrix(matrix, denominator))
        receipts.append(path)
    fresh = physical(x)
    matrices = []
    discrepancy = 0.0
    for key, a, ds, ph, q, ts, cols, *_ in blocks:
        scale = np.sqrt([sum(v*v for v in column.values()) for column in cols])
        old = np.array([projected_matrix(cols, action) for action in old_actions], float) / scale[:, None] / scale[None, :]
        diagonal = np.array([float(local_value(next(iter(column)), terms)) for column in cols])
        base = a + np.einsum('i,ijk->jk', np.array(coefficients, float), old) + np.diag(diagonal) + float(alpha)*ph + float(beta)*q
        direct = np.array(fresh[key][1], float) / scale[:, None] / scale[None, :]
        direct += np.einsum('i,ijk->jk', np.array(coefficients, float), old) + np.diag(diagonal) + float(alpha)*ph + float(beta)*q
        discrepancy = max(discrepancy, float(np.max(abs(direct - base))))
        new = np.array([projected_matrix(cols, action) for action in new_actions], float) / scale[:, None] / scale[None, :]
        matrices.append((key, base, new))
    if discrepancy > 1e-9:
        raise ValueError('Fresh physical matrix reconstruction disagrees')
    calls, cached, cached_result = 0, None, None
    best = None

    def evaluate(point):
        nonlocal calls, cached, cached_result, best
        point = np.asarray(point)
        if cached is not None and np.array_equal(cached, point):
            return cached_result
        if calls >= 250:
            raise RuntimeError('250 full-spectrum evaluations exhausted')
        calls += 1
        values, gradients = [], []
        for key, base, new in matrices:
            eigenvalues, vectors = eigh(base + np.einsum('i,ijk->jk', point, new), subset_by_index=[0, 0])
            vector = vectors[:, 0]
            values.append(float(eigenvalues[0]))
            gradients.append([float(vector @ matrix @ vector) for matrix in new])
        values, gradients = np.array(values), np.array(gradients)
        minimum = float(min(values))
        if max(abs(point)) <= 2 and (best is None or minimum > best[0]):
            best = (minimum, point.copy())
        cached, cached_result = point.copy(), (values, gradients)
        return cached_result

    baseline = float(min(evaluate(np.zeros(2))[0]))
    direction = np.array([.6, -.8])
    probe = np.array([.0013, -.0027])
    values, gradients = evaluate(probe)
    step = 1e-7
    derivative = (evaluate(probe + step*direction)[0] - evaluate(probe - step*direction)[0])/(2*step)
    # Some inactive symmetry blocks have exact degeneracies. Check the active minimum only.
    active = int(np.argmin(values))
    derivative_error = float(abs(derivative[active] - gradients[active] @ direction))
    if derivative_error > 2e-7:
        raise ValueError('Active spectral derivative check failed')
    initial = np.r_[np.zeros(2), baseline - 1e-8]
    constraint = {'type': 'ineq', 'fun': lambda z: evaluate(z[:2])[0] - z[2],
                  'jac': lambda z: np.column_stack([evaluate(z[:2])[1], -np.ones(len(matrices))])}
    outcome = None
    try:
        outcome = minimize(lambda z: -z[2], initial, jac=lambda z: np.array([0., 0., -1.]),
                           method='SLSQP', bounds=[(-2., 2.), (-2., 2.), (None, None)],
                           constraints=constraint, options={'maxiter': 150, 'ftol': 1e-12})
    except RuntimeError as exc:
        message = str(exc)
    rounded = [F(round(float(v)*10**6), 10**6) for v in best[1]]
    # Independent final spectral evaluation does not consume or evade an optimization iteration.
    final_minimum = min(float(eigh(base + np.einsum('i,ijk->jk', np.array(rounded, float), new),
                                   subset_by_index=[0, 0], eigvals_only=True)[0]) for key, base, new in matrices)
    ell = F(math.floor(final_minimum*10**7)-1, 10**7)
    lower = (ell - penalty_cost)/5
    accepted_lower = (F(certificate['penalized_lower']) - penalty_cost)/5
    proposal = {'accepted': False, 'kind': 'residual_coherence_fixed_recipe_proposal_v1',
                'seed_certificate': str(args.seed.resolve().relative_to(ROOT)),
                'new_coefficients': list(map(str, rounded)), 'proposed_penalized_lower': str(ell),
                'proposed_periodic_lower': str(lower), 'proposed_periodic_lower_float': float(lower),
                'proposed_gain_over_accepted_seed': str(lower - accepted_lower),
                'proposed_gain_over_accepted_seed_float': float(lower - accepted_lower),
                'numerical_baseline_periodic_lower': (baseline-float(penalty_cost))/5,
                'gain_over_same_numeric_baseline': (final_minimum-baseline)/5,
                'matrix_evaluations': calls + 1, 'max_optimizer_evaluations': 250,
                'optimizer_success': bool(outcome.success) if outcome is not None else False,
                'optimizer_message': outcome.message if outcome is not None else message,
                'derivative_error': derivative_error, 'fresh_physical_discrepancy': discrepancy,
                'local_blocks': len(matrices), 'seconds': time.monotonic()-started,
                'scope': 'Untrusted two-coordinate fixed-recipe spectral probe with coefficients bounded by2. All prior fields remain fixed. No production schema accepts this proposal and no exact PSD proof, family ceiling, convergence or transfer claim is made.'}
    files = {Path(__file__).resolve(), args.seed.resolve(), candidates_path, *receipts}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.is_relative_to(ROOT / 'experiments') or path.is_relative_to(BASE.parent / 'discovery'):
                files.add(path)
    proposal['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    (args.output / 'proposal.json').write_text(json.dumps(proposal, indent=2) + '\n')
    print(json.dumps({k: v for k, v in proposal.items() if k not in ('source_sha256', 'scope')}))


if __name__ == '__main__':
    main()
