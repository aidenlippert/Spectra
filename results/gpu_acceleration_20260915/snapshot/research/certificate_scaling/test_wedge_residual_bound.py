"""Independent fixed-N controls for wedge_residual_bound."""
from fractions import Fraction as F
from itertools import combinations
import unittest
import numpy as np

from research.certificate_scaling.wedge_residual_bound import residual_bounds


def act_word(word, state, m):
    amp = 1
    occ = set(i for i in range(m) if state >> i & 1)
    # operator product acts rightmost first
    for create, i in reversed(word):
        sign = (-1) ** sum(j < i for j in occ)
        if create:
            if i in occ: return None, 0
            occ.add(i)
        else:
            if i not in occ: return None, 0
            occ.remove(i)
        amp *= sign
    return sum(1 << i for i in occ), amp


def sector_matrix(residual, m, n):
    states = [sum(1 << i for i in c) for c in combinations(range(m), n)]
    idx = {s: i for i, s in enumerate(states)}
    out = np.zeros((len(states), len(states)))
    for word, coeff in residual.items():
        for j, s in enumerate(states):
            t, a = act_word(word, s, m)
            if t is not None: out[idx[t], j] += float(coeff) * a
    return out


class WedgeResidualControls(unittest.TestCase):
    def test_random_hermitian_bodies_bound_all_sectors(self):
        for k in (1, 2, 3):
            words = list(combinations(range(4), k))
            residual = {}
            for a, i in enumerate(words):
                for b, j in enumerate(words):
                    c = F((a + 2*b + 1) % 5 - 2)
                    if c:
                        # coefficient convention is the normal ordered body;
                        # enforce Hermiticity in the raw body matrix.
                        # residual_bounds applies the wedge reversal sign;
                        # choose raw coefficients so its coefficient matrix is
                        # the desired Hermitian matrix.
                        sign = (-1) ** (k * (k - 1) // 2)
                        residual[tuple((1, x) for x in i) + tuple((0, x) for x in j)] = c / sign
            # symmetrize exactly
            raw = residual
            residual = {}
            for i in words:
                for j in words:
                    w = tuple((1, x) for x in i) + tuple((0, x) for x in j)
                    wr = tuple((1, x) for x in j) + tuple((0, x) for x in i)
                    c = (raw.get(w, 0) + raw.get(wr, 0)) / 2
                    if c: residual[w] = c
            for n in range(5):
                lower, _ = residual_bounds(residual, 4, n)
                ev = np.linalg.eigvalsh(sector_matrix(residual, 4, n)).min()
                self.assertLessEqual(float(lower), ev + 1e-8)

    def test_positive_projector_has_nonnegative_bound(self):
        # (a_0+a_1)^dagger (a_0+a_1) is positive on every sector.
        # Its normal-ordered one-body expansion is a small exact control.
        residual = {((1, 0), (0, 0)): F(1), ((1, 1), (0, 1)): F(1),
                    ((1, 0), (0, 1)): F(1), ((1, 1), (0, 0)): F(1)}
        lower, _ = residual_bounds(residual, 2, 1)
        self.assertGreaterEqual(lower, 0)

    def test_nonhermitian_and_unbalanced_rejected(self):
        with self.assertRaises(ValueError):
            residual_bounds({((1, 0), (0, 1)): F(1)}, 3, 1)
        with self.assertRaises(ValueError):
            residual_bounds({((1, 0),): F(1)}, 3, 1)

    def test_k_greater_than_n_vanishes(self):
        word = tuple((1, i) for i in (0, 1)) + tuple((0, i) for i in (0, 1))
        total, details = residual_bounds({word: F(1)}, 4, 1)
        self.assertEqual(total, F(0))
        self.assertEqual(details[0]['fixed_N_lower'], '0')


if __name__ == '__main__':
    unittest.main()
