"""Exact orbital-coefficient lifts. No N-electron state-space construction.

Matrices are sparse dictionaries (increasing tuple I, increasing tuple J) ->
real rational. E_IJ has annihilators in reverse order; no factorial is used.
This implements and checks the supplied derivation, not an energy acceptance
shortcut. The unchanged CAR checker still re-expands every molecular proof.
"""
from fractions import Fraction as F


def checked(matrix, m, k):
    if type(m) is not int or type(k) is not int or not 0 <= k <= m:
        raise ValueError('Invalid orbital coefficient space')
    out = {}
    for key, value in matrix.items():
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError('Matrix keys must be two tuples')
        for indices in key:
            if (not isinstance(indices, tuple) or len(indices) != k or
                any(type(i) is not int or not 0 <= i < m for i in indices) or
                tuple(sorted(set(indices))) != indices):
                raise ValueError('Invalid increasing orbital tuple')
        if type(value) not in (int, F):
            raise ValueError('Exact real rational coefficients required')
        if value:
            out[key] = F(value)
    return out


def combine(*terms):
    out = {}
    for weight, matrix in terms:
        for key, value in matrix.items():
            out[key] = out.get(key, 0) + weight * value
    return {key: F(value) for key, value in out.items() if value}


def lift(matrix, m, k):
    matrix = checked(matrix, m, k)
    out = {}
    for (I, J), value in matrix.items():
        for r in range(m):
            if r in I or r in J:
                continue
            p = sum(i < r for i in I)
            q = sum(j < r for j in J)
            key = (I[:p] + (r,) + I[p:], J[:q] + (r,) + J[q:])
            out[key] = out.get(key, 0) + (-1)**(p+q) * value
    return {key: value for key, value in out.items() if value}


def contract(matrix, m, k):
    matrix = checked(matrix, m, k)
    if k == 0:
        raise ValueError('Cannot contract scalar')
    out = {}
    for (I, J), value in matrix.items():
        for p, r in enumerate(I):
            if r in J:
                q = J.index(r)
                key = (I[:p] + I[p+1:], J[:q] + J[q+1:])
                out[key] = out.get(key, 0) + (-1)**(p+q) * value
    return {key: value for key, value in out.items() if value}


def inverse_contracted(K, m):
    if type(m) is not int or m < 5:
        raise ValueError('Injective two-to-three lift requires m >= 5')
    K = checked(K, m, 2)
    T = contract(K, m, 2)
    t = sum(v for (I, J), v in T.items() if I == J) / F(3*(m-2))
    identity = {((i,), (i,)): F(1) for i in range(m)}
    B = combine((F(1, 2*(m-3)), T), (-t/F(2*(m-3)), identity))
    return combine((F(1, m-4), K), (-F(1, m-4), lift(B, m, 1)))


def inner(A, B):
    return sum((value * B.get(key, 0) for key, value in A.items()), F(0))


def split_three_body(W, m):
    W = checked(W, m, 3)
    K = contract(W, m, 3)
    V = inverse_contracted(K, m)
    Z = combine((1, W), (-1, lift(V, m, 2)))
    if contract(Z, m, 3):
        raise AssertionError('Contraction-free decomposition failed')
    if inner(Z, Z) != inner(W, W) - inner(V, K):
        raise AssertionError('Orthogonal norm identity failed')
    return V, Z


def factorized_lift_check(factors, m):
    """Signed triple-vector factors; builds K but never the triple matrix W."""
    parsed = []
    K = {}
    for weight, vector in factors:
        if type(weight) not in (int, F):
            raise ValueError('Exact signed weight required')
        clean = checked({(I, I): value for I, value in vector.items()}, m, 3)
        vector = {I: value for (I, J), value in clean.items()}
        parsed.append((F(weight), vector))
        for r in range(m):
            u = {I[:p]+I[p+1:]: (-1)**p * value
                 for I, value in vector.items() for p in range(3) if I[p] == r}
            for I, a in u.items():
                for J, b in u.items():
                    K[I, J] = K.get((I, J), 0) + weight*a*b
    K = {key: value for key, value in K.items() if value}
    norm = F(0)
    for t, (weight, vector) in enumerate(parsed):
        for u in range(t, len(parsed)):
            other_weight, other = parsed[u]
            dot = sum((v * other.get(I, 0) for I, v in vector.items()), F(0))
            norm += (1 if t == u else 2) * weight * other_weight * dot**2
    V = inverse_contracted(K, m)
    remainder = norm - inner(V, K)
    if norm < 0 or remainder < 0:
        raise AssertionError('Negative exact squared norm')
    return V, {'squared_norm': norm, 'remainder_squared_norm': remainder,
               'is_lift': remainder == 0, 'triple_matrix_materialized': False}


def from_polynomial(poly, m, k):
    """Canonical CAR annihilators ascend, giving (-1)^(k(k-1)/2)."""
    out = {}
    phase = (-1)**(k*(k-1)//2)
    for word, value in poly.items():
        if len(word) != 2*k:
            continue
        I = tuple(i for c, i in word if c)
        J = tuple(i for c, i in word if not c)
        if len(I) != k or len(J) != k or word != tuple((1, i) for i in I)+tuple((0, j) for j in J):
            raise ValueError('Expected canonical number-conserving CAR terms')
        out[I, J] = phase * value
    return checked(out, m, k)


def to_polynomial(matrix, m, k):
    phase = (-1)**(k*(k-1)//2)
    return {tuple((1, i) for i in I)+tuple((0, j) for j in J): phase*value
            for (I, J), value in checked(matrix, m, k).items()}
