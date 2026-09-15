"""Optional exact SPD witnesses for bounded integer PSD replay.

For invertible upper triangular integer R, strict positive diagonal dominance
of R^T A R proves A positive definite. Floating point never enters acceptance.
Matrices without witnesses retain the original singular-safe PSD elimination.
"""
import hashlib
import json

from experiments.marginal_polynomial_sos import integer_psd


def matrix_digest(matrix):
    return hashlib.sha256(json.dumps(matrix, separators=(',', ':')).encode()).hexdigest()


class CongruenceProofs:
    def __init__(self, witnesses):
        if type(witnesses) is not dict or len(witnesses) > 256:
            raise ValueError('At most 256 congruence witnesses required')
        for key, factor in witnesses.items():
            if type(key) is not str or len(key) != 64 or any(c not in '0123456789abcdef' for c in key):
                raise ValueError('Exact matrix SHA256 key required')
            if type(factor) is not list or not 1 <= len(factor) <= 300:
                raise ValueError('Bounded square integer congruence required')
            n = len(factor)
            if any(type(row) is not list or len(row) != n for row in factor):
                raise ValueError('Square integer congruence required')
            if any(type(v) is not int or abs(v) > 10**30 for row in factor for v in row):
                raise ValueError('Bounded integer congruence entries required')
            if any(factor[i][j] for i in range(n) for j in range(i)) or any(not factor[i][i] for i in range(n)):
                raise ValueError('Invertible upper triangular congruence required')
        self.witnesses = witnesses
        self.used = set()

    def check(self, matrix, initial_divisor=1):
        n = len(matrix)
        if not 1 <= n <= 300 or any(len(row) != n for row in matrix):
            raise ValueError('Bounded square moment matrix required')
        if any(type(v) is not int for row in matrix for v in row):
            raise ValueError('Integer moment matrix required')
        if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(n)):
            raise ValueError('Exact symmetric moment matrix required')
        if type(initial_divisor) is not int or initial_divisor <= 0:
            raise ValueError('Positive integer initial divisor required')
        key = matrix_digest(matrix)
        if key not in self.witnesses:
            return integer_psd(matrix, initial_divisor=initial_divisor)
        factor = self.witnesses[key]
        if len(factor) != n:
            raise ValueError('Congruence dimension differs from reconstructed matrix')
        product = [[sum(matrix[i][k]*factor[k][j] for k in range(j+1))
                    for j in range(n)] for i in range(n)]
        transformed = [[sum(factor[k][i]*product[k][j] for k in range(i+1))
                        for j in range(n)] for i in range(n)]
        margins = [row[i]-sum(abs(v) for j, v in enumerate(row) if i != j)
                   for i, row in enumerate(transformed)]
        if min(margins) <= 0:
            raise ValueError('Exact strict positive diagonal dominance failed')
        self.used.add(key)
        return {'dimension': n, 'rank': n, 'nullity': 0, 'positive_semidefinite': True,
                'method': 'integer_congruence_strict_diagonal_dominance',
                'matrix_sha256': key, 'minimum_exact_margin': str(min(margins))}

    def finish(self):
        if self.used != set(self.witnesses):
            raise ValueError('Unused congruence witnesses do not match reconstructed matrices')
        return {'supplied': len(self.witnesses), 'used': len(self.used)}
