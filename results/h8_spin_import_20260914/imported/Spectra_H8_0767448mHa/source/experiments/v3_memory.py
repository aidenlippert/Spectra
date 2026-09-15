"""Finite retained-coordinate error bounds for a fully dissipative block model.

Declared full generator:
  L = [[K - diag(retained_rates), B], [-B.T, -diag(rates)]].
K is skew, all rates positive, and the initial discarded vector is zero.
The bounds concern retained Euclidean coordinates, not physical observables.
"""
from fractions import Fraction as F
from math import isqrt


def _sqrt_up(value, digits=12):
    scale = 10 ** digits
    numerator = value.numerator * scale * scale
    denominator = value.denominator
    root = isqrt(numerator // denominator)
    if root * root * denominator < numerator: root += 1
    return F(root, scale)


def _checked(K, rates, B, initial, retained_rates):
    K = tuple(tuple(F(x) for x in row) for row in K)
    rates = tuple(F(x) for x in rates)
    B = tuple(tuple(F(x) for x in row) for row in B)
    initial = tuple(F(x) for x in initial)
    n, q = len(K), len(rates)
    if not 1 <= n <= 32 or not 1 <= q <= 32 or any(len(row) != n for row in K):
        raise ValueError('bounded square retained block and nonempty discarded block required')
    if len(B) != n or any(len(row) != q for row in B) or len(initial) != n:
        raise ValueError('incompatible block dimensions')
    if any(K[i][j] != -K[j][i] for i in range(n) for j in range(n)):
        raise ValueError('K must be skew')
    if any(r <= 0 for r in rates): raise ValueError('positive discarded rates required')
    retained_rates = tuple(F(x) for x in retained_rates) if retained_rates is not None else (min(rates),) * n
    if len(retained_rates) != n or any(r <= 0 for r in retained_rates):
        raise ValueError('positive retained rates required')
    norm1 = max(sum(abs(B[i][j]) for i in range(n)) for j in range(q))
    norminf = max(sum(abs(x) for x in row) for row in B)
    norm_b = _sqrt_up(norm1 * norminf)
    norm_initial = _sqrt_up(sum(x * x for x in initial))
    gamma = min(rates + retained_rates)
    return K, rates, B, initial, retained_rates, norm_b, norm_initial, gamma


def certificate(K, rates, B, initial, times=(1,), retained_rates=None):
    K, rates, B, initial, retained_rates, b, norm, gamma = _checked(K, rates, B, initial, retained_rates)
    times = tuple(F(t) for t in times)
    if not 1 <= len(times) <= 128 or any(t < 0 for t in times): raise ValueError('bounded nonnegative time list required')
    return {'K': [[str(x) for x in row] for row in K], 'rates': [str(x) for x in rates],
            'retained_rates': [str(x) for x in retained_rates],
            'B': [[str(x) for x in row] for row in B], 'initial': [str(x) for x in initial],
            'discarded_initial_is_zero': True,
            'B_norm_upper': str(b), 'gamma': str(gamma), 'initial_norm_upper': str(norm),
            'times': [{'time': str(t), 'finite_error_bound': str(b * b * t * t * norm / 2),
                       'uniform_error_bound': str(2 * b * b * norm / (gamma * gamma))} for t in times],
            'scope': 'retained coefficient 2-norm; full model as specified in module docstring'}


def verify(c):
    try:
        # Recheck model obligations without calling the certificate generator.
        K, rates, B, initial, retained, b, norm, gamma = _checked(c['K'], c['rates'], c['B'], c['initial'], c['retained_rates'])
        if c['discarded_initial_is_zero'] is not True: return False
        if F(c['B_norm_upper']) < b or F(c['initial_norm_upper']) < norm or F(c['gamma']) != gamma:
            return False
        if not isinstance(c['times'], list) or not 1 <= len(c['times']) <= 128: return False
        for row in c['times']:
            t = F(row['time'])
            if t < 0 or F(row['finite_error_bound']) < b * b * t * t * norm / 2:
                return False
            if F(row['uniform_error_bound']) < 2 * b * b * norm / (gamma * gamma):
                return False
        return True
    except (AttributeError, KeyError, TypeError, ValueError, IndexError, ZeroDivisionError):
        return False
