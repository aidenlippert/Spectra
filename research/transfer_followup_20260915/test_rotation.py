from fractions import Fraction as F
from itertools import combinations
import unittest
from experiments.marginal_symbolic import canonical, encode
from research.transfer_solver_20260915.independent_oracle import matrix
from research.transfer_followup_20260915.rotated_upper import transformed, rotate_and_round, mm, transpose, check_rotation, rounded_division


class Rotation(unittest.TestCase):
    def test_two_particle_unitary_with_independent_bit_oracle(self):
        Z, den = [[3, -4], [4, 3]], 5
        terms = {(): F(1, 7)}
        for i in range(4):
            terms[((1, i), (0, i))] = F(i+1, 3)
        terms[((1, 0), (0, 2))] = terms[((1, 2), (0, 0))] = F(1, 3)
        terms[((1, 0), (1, 1), (0, 3), (0, 2))] = F(2, 7)
        terms[((1, 2), (1, 3), (0, 1), (0, 0))] = F(2, 7)
        terms[((1, 0), (1, 2), (0, 2), (0, 0))] = F(5, 7)
        terms[((1, 1), (1, 3), (0, 3), (0, 1))] = F(2, 11)
        data = {'modes': 4, 'particles': 2, 'hamiltonian': encode(canonical(terms))}
        constant, one, two, hd, pairs = transformed(data, Z, den)
        exact = {(): F(constant, hd)}
        for i in range(4):
            for j in range(4):
                exact[((1, i), (0, j))] = F(one[i][j], hd*den**2)
        for a, (p, q) in enumerate(pairs):
            for b, (r, s) in enumerate(pairs):
                exact[((1, p), (1, q), (0, s), (0, r))] = F(two[a][b], hd*den**4)
        rotated = {**data, 'hamiltonian': encode(canonical(exact))}
        states, H = matrix(data)
        rotated_states, K = matrix(rotated)
        self.assertEqual(states, rotated_states)
        V = [[F(Z[i//2][j//2], den) if i % 2 == j % 2 else F(0) for j in range(4)] for i in range(4)]
        compound = [[V[p][i]*V[q][j]-V[p][j]*V[q][i] for i, j in pairs] for p, q in pairs]
        self.assertEqual(K, mm(mm(transpose(compound), H), compound))
        rounded, error = rotate_and_round(data, Z, den)
        _, R = matrix(rounded)
        self.assertTrue(all(sum(abs(a-b) for a, b in zip(row, other)) <= error for row, other in zip(K, R)))

    def test_nonorthogonal_rotation_refused(self):
        with self.assertRaisesRegex(ValueError, 'not exactly orthogonal'):
            check_rotation([[3, -4], [4, 2]], 5)

    def test_signed_ties_to_even(self):
        for n in range(-9, 10):
            self.assertEqual(rounded_division(n, 2), round(F(n, 2)))


if __name__ == '__main__':
    unittest.main()
