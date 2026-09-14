"""Sparse physical Krylov proposals with exact variational and residual replay."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_general_schur import gram
from experiments.marginal_sparse_response import q_action
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_schur_transfer import ldl_pivots


def witness_statistics(oracle, witness):
    """Recompute the normalized Rayleigh value and full eigen-residual exactly."""
    upper = oracle.upper(witness)
    vector = {s: F(x) for s, x in zip(witness['states'], witness['amplitudes']) if x}
    image = q_action(oracle, set(), vector)
    norm = sum(x*x for x in vector.values())
    residual = {s: image.get(s, F(0)) - upper * vector.get(s, F(0)) for s in set(image) | set(vector)}
    variance = sum(x*x for x in residual.values()) / norm
    return upper, variance, len(vector)


def replay(certificate):
    if certificate.get('kind') != 'sparse_variational_upper_v1':
        raise ValueError('Unsupported sparse upper certificate')
    oracle = DeterminantOracle(certificate)
    upper, variance, support = witness_statistics(oracle, certificate.get('independent_upper'))
    return {'upper': rational_text(upper), 'upper_float': float(upper),
            'residual_norm_squared': rational_text(variance), 'residual_norm_squared_float': float(variance),
            'support_size': support, 'unique_action_states': len(oracle.cache),
            'referenced_determinants': oracle.referenced_state_count(),
            'scope': 'Exact variational upper and eigen-residual for the rational Hamiltonian; no global lower endpoint or ground-state error bound'}


def refine(fixture_path, proposal_path, out, steps=12):
    import numpy as np
    if type(steps) is not int or not 1 <= steps <= 32:
        raise ValueError('Physical Krylov budget must be one to 32')
    out = Path(out)
    if out.exists():
        raise ValueError('Preserve previous sparse upper export')
    fixture = json.loads(Path(fixture_path).read_text())
    proposal = json.loads(Path(proposal_path).read_text())
    oracle = DeterminantOracle(fixture)
    witness = proposal['independent_upper']
    initial = oracle.upper(witness)
    first = {s: F(x) for s, x in zip(witness['states'], witness['amplitudes']) if x}
    basis, images, history = [first], [], []
    best = dict(witness)
    best_upper = initial
    started = time.monotonic()
    for iteration in range(steps):
        images.append(q_action(oracle, set(), basis[-1]))
        metric = gram(basis, basis)
        if ldl_pivots(metric) is None:
            raise ValueError('Physical Krylov proposal basis is dependent')
        projected = gram(basis, images)
        chol = np.linalg.cholesky(np.array(metric, dtype=float))
        left = np.linalg.solve(chol, np.array(projected, dtype=float))
        normalized = np.linalg.solve(chol, left.T).T
        _, vectors = np.linalg.eigh(normalized)
        coefficients = np.linalg.solve(chol.T, vectors[:, 0])
        physical = {}
        for column, coefficient in zip(basis, coefficients):
            for s, x in column.items():
                physical[s] = physical.get(s, 0.) + float(x) * coefficient
        integers = {s: int(round(x * 10**10)) for s, x in physical.items()}
        states = sorted(s for s, x in integers.items() if x)
        candidate = {'states': states, 'amplitudes': [integers[s] for s in states]}
        energy = oracle.upper(candidate)
        if energy < best_upper:
            best, best_upper = candidate, energy
        history.append({'basis_dimension': len(basis), 'upper': rational_text(best_upper),
                        'support_size': sum(bool(x) for x in best['amplitudes'])})
        if iteration + 1 == steps:
            break
        next_column = {s: float(x) for s, x in images[-1].items()}
        for _ in range(2):
            for column in basis:
                norm = sum(float(x)**2 for x in column.values())
                overlap = sum(float(x) * next_column.get(s, 0.) for s, x in column.items()) / norm
                for s, x in column.items():
                    next_column[s] = next_column.get(s, 0.) - overlap * float(x)
        norm = math.sqrt(sum(x*x for x in next_column.values()))
        if norm <= 1e-12:
            break
        integers = {s: int(round(x / norm * 10**10)) for s, x in next_column.items()}
        basis.append({s: F(x) for s, x in integers.items() if x})
    c = {key: fixture[key] for key in ('modes', 'particles', 'hamiltonian')}
    c.update(kind='sparse_variational_upper_v1', independent_upper=best)
    receipt = replay(c)
    receipt.update(initial_upper=rational_text(initial), improvement=rational_text(initial - best_upper),
                   basis_dimension=len(images), elapsed_seconds=time.monotonic() - started,
                   construction_unique_action_states=len(oracle.cache),
                   construction_referenced_determinants=oracle.referenced_state_count(),
                   source_proposal=str(proposal_path))
    out.mkdir(parents=True)
    for filename, value in (('certificate.json', c), ('receipt.json', receipt), ('history.json', history)):
        (out / filename).write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ('upper_float', 'residual_norm_squared_float', 'support_size', 'basis_dimension', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture')
    parser.add_argument('--proposal')
    parser.add_argument('--output')
    parser.add_argument('--steps', type=int, default=12)
    parser.add_argument('--verify')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.fixture and args.proposal and args.output:
        refine(args.fixture, args.proposal, args.output, args.steps)
    else:
        parser.error('Supply --fixture, --proposal and --output, or --verify')
