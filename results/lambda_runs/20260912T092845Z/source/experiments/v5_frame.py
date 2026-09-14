"""Exact dense mode geometry for the bounded scalar-sensing construction."""
from fractions import Fraction as F

# Columns of a real orthogonal Hadamard transformation; laboratory coordinates
# are dense mixtures of the coordinates used in the proof.
Q = tuple(tuple(F(value, 2) for value in column) for column in
          ((1, 1, 1, 1), (1, -1, 1, -1), (1, 1, -1, -1), (1, -1, -1, 1)))
BASE_LAMBDAS = (10**14, 10**10, 10**6)
C, S = F(3, 5), F(4, 5)


def dot(a, b):
    if len(a) != 4 or len(b) != 4:
        raise ValueError('four-dimensional vectors required')
    return sum((x * y for x, y in zip(a, b)), F(0))


def _combine(a, x, b, y):
    return tuple(a * xi + b * yi for xi, yi in zip(x, y))


def _check(signs):
    if type(signs) is not tuple or len(signs) > 3 or any(type(s) is not int or s not in (-1, 1) for s in signs):
        raise ValueError('strict sign tuple of length 0..3 required')


def frame(signs):
    _check(signs)
    residual, modes = Q[0], []
    for index, sign in enumerate(signs, 1):
        new_axis = Q[index]
        mode = _combine(C, residual, sign * S, new_axis)
        residual = _combine(-sign * S, residual, C, new_axis)
        modes.append(mode)
    return tuple(modes), residual


def null_probe(known_prefix):
    """Null all known modes and the next +1 hypothesis."""
    _check(known_prefix)
    if len(known_prefix) == 3:
        raise ValueError('no fourth discovery in this family')
    _, residual = frame(known_prefix)
    return _combine(-S, residual, C, Q[len(known_prefix) + 1])


def covariance(signs, alpha=F(1)):
    _check(signs)
    if not signs:
        raise ValueError('at least one active mode required')
    if type(alpha) is not F or not 1 <= alpha <= 2:
        raise ValueError('exact coupling scale in [1,2] required')
    modes, _ = frame(signs)
    return tuple(tuple(F(int(i == j)) + sum(alpha * lam * mode[i] * mode[j]
                 for mode, lam in zip(modes, BASE_LAMBDAS))
                 for j in range(4)) for i in range(4))
