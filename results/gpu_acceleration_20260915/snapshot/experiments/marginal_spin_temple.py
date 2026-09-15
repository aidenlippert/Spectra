"""Exact spin-sector eigenvalue counting followed by a Temple energy interval.

Only the lowest energy transfers through SU(2), not the excited-state count.
Numerical prefix selection is optional discovery; replay uses rational arithmetic.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_sparse_response import parse_basis, prepare, response_matrix, propose_direction
from experiments.marginal_sparse_upper import witness_statistics
from experiments.marginal_spin_reduction import SpinZeroOracle, spin_gap
from experiments.marginal_symbolic import add, scale


def inertia(matrix):
    """Exact congruence with 1x1/2x2 pivots; zero eigenvalues are counted."""
    n = len(matrix)
    if not 1 <= n <= 32 or any(len(row) != n for row in matrix):
        raise ValueError('Square inertia matrix of dimension one to 32 required')
    a = [[F(x) for x in row] for row in matrix]
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Exact symmetric inertia matrix required')
    result = {'positive': 0, 'negative': 0, 'zero': 0}
    while a:
        n = len(a)
        pivot = next((i for i in range(n) if a[i][i]), None)
        if pivot is not None:
            d = a[pivot][pivot]
            result['positive' if d > 0 else 'negative'] += 1
            rest = [i for i in range(n) if i != pivot]
            a = [[a[i][j] - a[i][pivot]*a[pivot][j]/d for j in rest] for i in rest]
            continue
        pair = next(((i, j) for i in range(n) for j in range(i+1, n) if a[i][j]), None)
        if pair is None:
            result['zero'] += n
            break
        # All diagonals vanish: this nonsingular 2x2 block has inertia (1,1,0).
        p, q = pair
        d = a[p][q]
        rest = [i for i in range(n) if i not in pair]
        a = [[a[i][j] - (a[i][p]*a[q][j] + a[i][q]*a[p][j])/d
              for j in rest] for i in rest]
        result['positive'] += 1
        result['negative'] += 1
    return result


def replay(certificate):
    if certificate.get('kind') != 'spin_temple_interval_v1':
        raise ValueError('Unsupported spin Temple certificate')
    original = DeterminantOracle(certificate)
    reduced = certificate.get('spin_symmetric_certificate')
    if (type(reduced) is not dict or reduced.get('modes') != original.modes
            or reduced.get('particles') != original.particles):
        raise ValueError('Spin reference must describe the same physical sector')
    oracle = SpinZeroOracle(reduced)
    p = reduced.get('retained_states')
    retained = oracle.retained(p)
    gamma, tau = F(reduced['complement_lower']), F(reduced['excitation_lower'])
    if gamma <= tau:
        raise ValueError('Strict complement positivity at the excitation threshold required')
    coverage, reference = spin_gap(oracle, p, reduced)
    recipe = reduced.get('response_basis')
    basis = [] if recipe == [] else parse_basis(oracle, retained, recipe, 32)
    counted = inertia(response_matrix(prepare(oracle, p, gamma, basis), tau))
    if counted['negative'] > 1:
        raise ValueError('Schur lower matrix does not certify at most one eigenvalue below threshold')
    witness = reduced.get('independent_upper')
    mu, variance, support = witness_statistics(oracle, witness)
    if mu >= tau:
        raise ValueError('Exact witness Rayleigh value must lie below excitation threshold')
    # The Rayleigh witness forces one eigenvalue below tau. A zero count would
    # contradict the checked Schur lower bound and is never accepted silently.
    if counted['negative'] != 1:
        raise ValueError('Inconsistent spectral count and physical witness')
    symmetric_lower = mu - variance/(tau-mu)
    error = sum(abs(x) for x in add(original.h, scale(oracle.h, -1)).values())
    lower, upper = symmetric_lower - error, original.upper(witness)
    if lower > upper:
        raise ValueError('Inconsistent original-Hamiltonian Temple interval')
    oracles = [original, oracle] + ([reference] if reference is not None else [])
    sources = set().union(*(set(o.cache) for o in oracles))
    referenced = sources | {s for o in oracles for row in o.cache.values() for s in row}
    result = {'lower': rational_text(lower), 'upper': rational_text(upper),
        'width': rational_text(upper-lower), 'width_float': float(upper-lower),
        'spin_symmetric_lower': rational_text(symmetric_lower),
        'spin_symmetric_upper': rational_text(mu),
        'spin_symmetric_width': rational_text(mu-symmetric_lower),
        'spin_symmetric_width_float': float(mu-symmetric_lower),
        'spin_sector_excitation_lower': rational_text(tau),
        'excitation_separation_from_witness': rational_text(tau-mu),
        'residual_norm_squared': rational_text(variance), 'residual_norm_squared_float': float(variance),
        'spin_symmetry_error_bound': rational_text(error), 'schur_lower_inertia': counted,
        'complement_coverage': coverage, 'retained_dimension': len(p),
        'response_dimension': len(basis), 'witness_support': support,
        'spin_sector_dimension': oracle.sector_dimension, 'full_sector_dimension': original.sector_dimension,
        'unique_determinant_sources': len(sources), 'referenced_determinants': len(referenced),
        'spin_symmetric_action_states': len(oracle.cache), 'original_action_states': len(original.cache),
        'scope': 'Exact full-P Schur lower inertia and strict Q gap certify the second eigenvalue only in Sz=0 of Hs. Temple uses the recomputed Rayleigh value and full variance. SU(2) transfers the ground minimum only; original-H norm transfer is separate. No full-sector excited gap, configuration compression, or scaling claim.'}
    if reference is not None:
        result['reference_action_states'] = len(reference.cache)
    return result


def prefix_data(data, count):
    result = dict(data)
    for key in ('metric', 'hamiltonian', 'squared'):
        result[key] = [row[:count] for row in data[key][:count]]
    for key in ('coupling', 'h_coupling'):
        result[key] = data[key][:count]
    return result


def excitation_threshold(mu, gamma, offset):
    offset = F(offset)
    if offset <= 0:
        raise ValueError('Positive excitation offset required')
    tau = F((mu+offset)*10**12 // 1, 10**12)
    if not mu < tau < gamma:
        raise ValueError('Excitation proposal must lie strictly between witness and complement bound')
    return tau


def discover(original, reduced, oracle, offset, target):
    """Fresh response proposals at one excitation threshold, with exact final replay.

    Selecting the second Schur eigenvector is a bounded heuristic, not a
    convergence theorem. Failure at the response cap returns no certificate.
    """
    import numpy as np
    p = reduced['retained_states']
    retained = oracle.retained(p)
    gamma = F(reduced['complement_lower'])
    mu, _, _ = witness_statistics(oracle, reduced['independent_upper'])
    tau = excitation_threshold(mu, gamma, offset)
    recipe, basis, history = [], [], []
    for count in range(33):
        data = prepare(oracle, p, gamma, basis)
        values, vectors = np.linalg.eigh(response_matrix(data, float(tau), True))
        negatives = int(np.count_nonzero(values < 0))
        history.append({'response_dimension': count, 'numerical_negative_count': negatives,
                        'smallest_eigenvalues': [float(x) for x in values[:2]]})
        if negatives <= 1:
            exact = inertia(response_matrix(data, tau))
            history[-1]['exact_inertia'] = exact
            if exact['negative'] == 1:
                candidate = dict(reduced, response_basis=recipe, excitation_lower=rational_text(tau))
                candidate.pop('kind', None)
                candidate.pop('lower', None)
                final = dict(original, kind='spin_temple_interval_v1', spin_symmetric_certificate=candidate)
                receipt = replay(final)
                if F(receipt['width']) > target:
                    raise ValueError('Witness variance and excitation threshold do not meet the energy target')
                return final, receipt, history
        if count == 32:
            break
        direction = 1 if len(p) > 1 else 0
        column, steps = propose_direction(oracle, retained, data['w'], vectors[:, direction], float(tau), basis)
        recipe.append(column)
        basis = parse_basis(oracle, retained, recipe, 32)
        history[-1]['next_proposal_cg_steps'] = steps
    raise ValueError('Fresh excitation response exhausted the 32-direction budget')


def compress(source, output, offset=F(1, 1000)):
    """Reuse proposed response prefixes; independently verify all final premises."""
    import numpy as np
    source, out, offset = Path(source), Path(output), F(offset)
    if out.exists():
        raise ValueError('Preserve previous Temple export')
    old = json.loads(source.read_text())
    if old.get('kind') != 'spin_reduced_perturbation_interval_v1':
        raise ValueError('Expected a spin-reduced response source')
    started = time.monotonic()
    DeterminantOracle(old)
    reduced = dict(old['spin_symmetric_certificate'])
    oracle = SpinZeroOracle(reduced)
    p = reduced['retained_states']
    basis = parse_basis(oracle, oracle.retained(p), reduced['response_basis'], 32)
    mu, _, _ = witness_statistics(oracle, reduced['independent_upper'])
    # Short rational threshold; all strict comparisons are checked again exactly.
    gamma = F(reduced['complement_lower'])
    tau = excitation_threshold(mu, gamma, offset)
    data = prepare(oracle, p, gamma, basis)
    history = []
    for count in range(len(basis)+1):
        values = np.linalg.eigvalsh(response_matrix(prefix_data(data, count), float(tau), True))
        numerical_negative = int(np.count_nonzero(values < 0))
        history.append({'response_dimension': count, 'numerical_negative_count': numerical_negative,
                        'smallest_eigenvalues': [float(x) for x in values[:2]]})
        if numerical_negative > 1:
            continue
        exact = inertia(response_matrix(prefix_data(data, count), tau))
        history[-1]['exact_inertia'] = exact
        if exact['negative'] == 1:
            break
    else:
        raise ValueError('No inherited prefix certifies the excitation count')
    reduced['response_basis'] = reduced['response_basis'][:count]
    reduced['excitation_lower'] = rational_text(tau)
    reduced.pop('lower', None)
    reduced.pop('kind', None)
    certificate = {key: old[key] for key in ('modes', 'particles', 'hamiltonian')}
    certificate.update(kind='spin_temple_interval_v1', spin_symmetric_certificate=reduced)
    receipt = replay(certificate)
    receipt.update(source=str(source), elapsed_seconds=time.monotonic()-started,
        inherited_response_dimension=len(basis), construction_spin_symmetric_action_states=len(oracle.cache),
        discovery_scope='Prefix selection inherits a previously constructed response basis and witness; old energy endpoints are unused. Fewer final directions does not establish cheaper discovery.')
    out.mkdir(parents=True)
    for filename, value in (('certificate.json', certificate), ('receipt.json', receipt), ('prefix_scan.json', history)):
        (out/filename).write_text(json.dumps(value, indent=2)+'\n')
    print(json.dumps({key: receipt[key] for key in ('width_float', 'spin_symmetric_width_float', 'response_dimension', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify')
    parser.add_argument('--source')
    parser.add_argument('--output')
    parser.add_argument('--offset', default='1/1000')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.source and args.output:
        compress(args.source, args.output, F(args.offset))
    else:
        parser.error('Specify --verify or --source and --output')
