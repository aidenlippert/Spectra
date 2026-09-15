from fractions import Fraction as F
from itertools import combinations, permutations
import unittest
import numpy as np
from research.transfer_solver_20260915.independent_oracle import matrix
from research.transfer_followup_20260915.correlated_guide import compound, canonical_moment


def expectation(word, vector):
    dagger = tuple((1-a, i) for a, i in reversed(word))
    data = {'modes': 6, 'particles': 3,
            'hamiltonian': [{'word': word, 'coefficient': '1/2'}, {'word': dagger, 'coefficient': '1/2'}]}
    states, H = matrix(data)
    A = np.array(H, dtype=float)
    return float(vector@A@vector/(vector@vector))


class Transport(unittest.TestCase):
    def test_all_three_orders_against_independent_fock_action(self):
        Z = [[3, -4, 0], [4, 3, 0], [0, 0, 5]]
        Oexact = [[F(Z[i//2][j//2], 5) if i % 2 == j % 2 else F(0) for j in range(6)] for i in range(6)]
        O = np.array(Oexact, dtype=float)
        sets3 = list(combinations(range(6), 3))
        v = np.array([i-8 if sum(j % 2 == 0 for j in x) == 1 else 0 for i, x in enumerate(sets3)], dtype=float)
        # Independent explicit Leibniz determinant expansion for the state.
        W = []
        for a in sets3:
            row = []
            for b in sets3:
                total = F(0)
                for perm in permutations(range(3)):
                    sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
                    value = F(sign)
                    for i in range(3):
                        value *= Oexact[a[i]][b[perm[i]]]
                    total += value
                row.append(float(total))
            W.append(row)
        rotated = np.array(W)@v
        for degree in range(1, 4):
            for alpha in range(degree+1):
                sets = [x for x in combinations(range(6), degree) if sum(i % 2 == 0 for i in x) == alpha]
                if not sets:
                    continue
                words = [[tuple((1, i) for i in a)+tuple((0, i) for i in reversed(b)) for b in sets] for a in sets]
                D = np.array([[expectation(w, v) for w in row] for row in words])
                actual = compound(O, sets)@D@compound(O, sets).T
                expected = np.array([[expectation(w, rotated) for w in row] for row in words])
                np.testing.assert_allclose(actual, expected, atol=1e-13)
                key = degree, alpha
                lookup = {x: i for i, x in enumerate(sets)}
                w = tuple((1, i) for i in sets[0])+tuple((0, i) for i in sets[-1])
                self.assertAlmostEqual(canonical_moment(w, {key: actual}, {key: lookup}), expectation(w, rotated), places=12)


if __name__ == '__main__':
    unittest.main()
