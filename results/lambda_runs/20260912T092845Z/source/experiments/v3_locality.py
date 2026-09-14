"""Uniform physical-observable certificates for a restricted dissipative family.

No full Pauli space is constructed to generate or verify the certificate.
The model assumes one- and two-site Hamiltonian terms and known local
Pauli depolarization at rate gamma per nonidentity factor.
"""
from fractions import Fraction as F
from math import ceil
from pathlib import Path
from itertools import product
import json
import numpy as np
from experiments.pauli import commutator_i


def validate(hamiltonian, seed, gamma):
    if not isinstance(seed, str) or not 1 <= len(seed) <= 64 or any(p not in 'IXYZ' for p in seed) or sum(p != 'I' for p in seed) != 1:
        raise ValueError('one-site Pauli seed on 1..64 qubits required')
    if not isinstance(gamma, F) or gamma <= 0:
        raise ValueError('positive exact depolarization rate required')
    n = len(seed)
    incident = [F(0)] * n
    for label, value in hamiltonian.items():
        if not isinstance(label, str) or len(label) != n or any(p not in 'IXYZ' for p in label):
            raise ValueError('invalid Hamiltonian Pauli term')
        if not isinstance(value, F):
            raise ValueError('exact rational Hamiltonian required')
        support = [i for i, p in enumerate(label) if p != 'I']
        if not 1 <= len(support) <= 2:
            raise ValueError('the certificate admits only one- and two-site terms')
        for i in support:
            incident[i] += abs(value)
    return max(incident)


def neighborhood(hamiltonian, seed, depth, cap=512):
    if type(depth) is not int or not 0 <= depth <= 12 or type(cap) is not int or not 1 <= cap <= 4096:
        raise ValueError('finite depth/cap required')
    basis, frontier = {seed}, {seed}
    for _ in range(depth):
        new = set()
        for word in frontier:
            new.update(commutator_i(hamiltonian, word))
        new -= basis
        if len(basis) + len(new) > cap:
            return tuple(sorted(basis)), False
        basis.update(new)
        frontier = new
    return tuple(sorted(basis)), True


def exp_negative_interval(x, tolerance=F(1, 10**14)):
    """Exact enclosure via positive exp series and a geometric tail bound."""
    if not isinstance(x, F) or not 0 <= x <= 128 or tolerance <= 0:
        raise ValueError('bounded nonnegative exponent required')
    total, term, n = F(1), F(1), 0
    while n < 1024:
        ratio = x / (n + 2)
        next_term = term * x / (n + 1)
        if ratio < 1:
            tail = next_term / (1 - ratio)
            low, high = 1 / (total + tail), 1 / total
            if high - low <= tolerance:
                return low, high, n
        n += 1
        term = next_term
        total += term
    raise ValueError('exponential enclosure resource cap exceeded')


def certify(hamiltonian, seed, gamma=F(300), tolerance=F(1, 1000), cutoff=F(1, 32), max_depth=6, cap=512):
    kappa = validate(hamiltonian, seed, gamma)
    if not isinstance(cutoff, F) or cutoff <= 0 or not isinstance(tolerance, F) or tolerance <= 0:
        raise ValueError('positive exact tolerance and cutoff required')
    if type(max_depth) is not int or not 0 <= max_depth <= 12:
        raise ValueError('bounded maximum depth required')
    lam = gamma - 2 * kappa
    x = 2 * kappa * cutoff
    if lam <= 0 or x >= 1:
        return {'status': 'no_certificate', 'reason': 'damping or short-time condition fails'}
    _, exp_upper, terms = exp_negative_interval(lam * cutoff)
    late = 2 * exp_upper
    for depth in range(max_depth + 1):
        early = x ** (depth + 1) / (1 - x)
        bound = max(early, late)
        if bound <= tolerance:
            basis, complete = neighborhood(hamiltonian, seed, depth, cap)
            if not complete:
                return {'status': 'no_certificate', 'reason': 'representation cap reached'}
            return {'status': 'certified', 'seed': seed, 'basis': list(basis), 'depth': depth,
                    'gamma': str(gamma), 'kappa': str(kappa), 'lambda': str(lam),
                    'cutoff': str(cutoff), 'tolerance': str(tolerance),
                    'early_error_bound': str(early), 'late_error_bound': str(late),
                    'uniform_all_time_error_bound': str(bound), 'exp_series_terms': terms,
                    'scope': 'all density states; fixed admitted Hamiltonian; exact projected semigroup'}
    return {'status': 'no_certificate', 'reason': 'depth budget insufficient for requested error'}


def verify(hamiltonian, c):
    """Independent certificate obligations; does not call certify."""
    try:
        if c['status'] != 'certified': return False
        seed, gamma, cutoff = c['seed'], F(c['gamma']), F(c['cutoff'])
        kappa = validate(hamiltonian, seed, gamma)
        lam = gamma - 2 * kappa
        depth = c['depth']
        if type(depth) is not int or not 0 <= depth <= 12 or not 0 < cutoff <= 1 or lam <= 0:
            return False
        basis = c['basis']
        if not isinstance(basis, list) or not 1 <= len(basis) <= 4096 or len(set(basis)) != len(basis): return False
        if seed not in basis or any(len(w) != len(seed) or any(p not in 'IXYZ' for p in w) or set(w) == {'I'} for w in basis): return False
        required, complete = neighborhood(hamiltonian, seed, depth, 4096)
        if not complete or not set(required) <= set(basis): return False
        x = 2 * kappa * cutoff
        if x >= 1: return False
        early = x ** (depth + 1) / (1 - x)
        _, exp_upper, terms = exp_negative_interval(lam * cutoff)
        late = 2 * exp_upper
        fields = {'kappa': kappa, 'lambda': lam, 'early_error_bound': early,
                  'late_error_bound': late, 'uniform_all_time_error_bound': max(early, late)}
        return (all(c[k] == str(v) for k, v in fields.items())
                and c['exp_series_terms'] == terms and 0 < F(c['tolerance'])
                and max(early, late) <= F(c['tolerance']))
    except (AttributeError, KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def reduced_columns(hamiltonian, basis, gamma):
    where = {w: i for i, w in enumerate(basis)}
    columns = []
    for i, word in enumerate(basis):
        column = {where[w]: a for w, a in commutator_i(hamiltonian, word).items() if w in where}
        column[i] = -gamma * sum(p != 'I' for p in word)
        columns.append(column)
    return columns


def propagate_enclosed(hamiltonian, c, time, numerical_tolerance=F(1, 10**10)):
    """Rational uniformization of the small projected coefficient generator.

    Returns a rational center and an l1 numerical error bound, independently
    composable with the physical truncation certificate.
    """
    if not verify(hamiltonian, c): raise ValueError('invalid structural certificate')
    if not isinstance(time, F) or time < 0: raise ValueError('exact nonnegative time required')
    basis, gamma = c['basis'], F(c['gamma'])
    columns = reduced_columns(hamiltonian, basis, gamma)
    rate = max(-col[j] for j, col in enumerate(columns))
    x = rate * time
    if x > 64: raise ValueError('propagation exponent budget exceeded')
    for j, column in enumerate(columns):
        # M=I+A/rate has induced l1 norm <=1 under the certified damping condition.
        if abs(1 + column[j] / rate) + sum(abs(v / rate) for i, v in column.items() if i != j) > 1:
            raise AssertionError('uniformization is not contractive')
    power = [F(int(w == c['seed'])) for w in basis]
    partial = power.copy()
    weight, scalar_sum, m = F(1), F(1), 0
    while m < 512:
        ratio = x / (m + 2)
        next_weight = weight * x / (m + 1)
        if ratio < 1:
            remainder = next_weight / (1 - ratio)
            low, high = 1 / (scalar_sum + remainder), 1 / scalar_sum
            # The uncomputed Poisson tail <= high*remainder; scalar uncertainty
            # times the retained vector is bounded by its exact l1 norm.
            error = (high - low) * sum(abs(v) for v in partial) / 2 + high * remainder
            if error <= numerical_tolerance:
                midpoint = (low + high) / 2
                return {'center': [v * midpoint for v in partial], 'numerical_error_bound': error,
                        'series_terms': m, 'basis': basis}
        updated = power.copy()
        for j, column in enumerate(columns):
            for i, value in column.items():
                updated[i] += value * power[j] / rate
        power = updated
        m += 1
        weight = next_weight
        scalar_sum += weight
        partial = [a + weight * b for a, b in zip(partial, power)]
    raise ValueError('uniformization term budget exceeded')


def chain_hamiltonian(n, coupling=F(1, 5)):
    if type(n) is not int or not 2 <= n <= 64: raise ValueError('n must be in 2..64')
    def word(sites): return ''.join(sites.get(i, 'I') for i in range(n))
    h = {word({i: 'Z'}): F(1, 2) for i in range(n)}
    for i in range(n - 1):
        h[word({i: 'X', i + 1: 'X'})] = F(1)
        h[word({i: 'Z', i + 1: 'Z'})] = coupling
    return h, word({n // 2: 'Z'})


def run():
    rows = []
    for n in (3, 6, 12, 24, 48):
        h, seed = chain_hamiltonian(n)
        c = certify(h, seed)
        if not verify(h, c): raise AssertionError('certificate failed')
        row = {'qubits': n, 'hamiltonian_terms': len(h), 'representation_size': len(c['basis']),
               'full_nonidentity_pauli_space': 4 ** n - 1, 'certificate': c,
               'bound_float': float(F(c['uniform_all_time_error_bound']))}
        # One exact short-time reduced prediction per size, including numerical error.
        pred = propagate_enclosed(h, c, F(1, 1000))
        row['rational_prediction'] = {'time': '1/1000', 'numerical_error_bound': str(pred['numerical_error_bound']),
                                       'series_terms': pred['series_terms'],
                                       'coefficients': [str(x) for x in pred['center']]}
        axes = {n // 2: 'Y', n // 2 - 1: 'X', n // 2 + 1: 'X'}
        expectation = sum(value for word, value in zip(c['basis'], pred['center'])
                          if all(p == 'I' or axes.get(i) == p for i, p in enumerate(word)))
        error = F(c['uniform_all_time_error_bound']) + pred['numerical_error_bound']
        row['damping_only_counterexample'] = {
            'product_state_axes': axes, 'other_sites': 'maximally mixed',
            'projected_expectation': str(expectation),
            'full_expectation_interval': [str(expectation - error), str(expectation + error)],
            'damping_only_prediction': '0',
            'absolute_error_lower_bound': str(max(F(0), expectation - error)),
            'exceeds_requested_tolerance': expectation - error > F(c['tolerance'])}
        rows.append(row)
    result = {'rows': rows, 'claim': 'conditional compact approximate dynamics under strong known depolarization',
              'compounding_learning_demonstrated': False, 'physical_experiments_performed': False,
              'important_limit': 'Hamiltonian and dissipation model supplied; basis-selection procedure programmed; no stateful-baseline advantage'}
    out = Path(__file__).resolve().parents[1] / 'results/v3'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'locality_results.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    print(json.dumps([{k: v for k, v in row.items() if k not in ('certificate', 'rational_prediction')}
                      for row in run()['rows']], indent=2))
