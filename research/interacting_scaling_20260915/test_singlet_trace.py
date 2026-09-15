from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import add, adj, canonical, mono, product, number_shift, encode
from research.collective_completion_20260914.spin_screen import spin_squared, ladder_ideal
from research.certificate_scaling.commutator_dual_witness import psd
from research.interacting_scaling_20260915.singlet_trace import singlet_dimension, singlet_trace, magnetic_trace


class SingletTraceTests(unittest.TestCase):
    def test_normalization_number_and_singlet_ideals(self):
        for m, n in [(4, 2), (8, 4), (16, 8), (20, 10)]:
            self.assertEqual(singlet_trace(mono(()), m, n), 1)
            self.assertEqual(singlet_trace(spin_squared(m), m, n), 0)
            self.assertEqual(singlet_trace(product(number_shift(m, n), mono(((1, 0), (0, 0)))), m, n), 0)
            self.assertEqual(singlet_trace(ladder_ideal(m, mono(((1, 1), (0, 2)))), m, n), 0)
        self.assertEqual(singlet_dimension(16, 8), 1764)

    def test_twirl_is_necessary_for_arbitrary_operators(self):
        # In the two-orbital two-electron singlet sector, triplet magnetic
        # contributions must cancel after averaging this non-invariant A.
        p = mono(((1, 0), (0, 0)))
        raw = (magnetic_trace(p, 4, 1, 1)-magnetic_trace(p, 4, 2, 0))/3
        self.assertEqual(raw, F(1, 3))
        self.assertEqual(singlet_trace(p, 4, 2), F(1, 2))

    def test_gram_psd_and_forced_spin_nullspace(self):
        m, n = 8, 4
        raising = add(*(mono(((1, i), (0, i+1))) for i in range(0, m, 2)))
        polys = [raising, mono(()), mono(((1, 0), (0, 2)))]
        gram = [[singlet_trace(product(adj(p), q), m, n) for q in polys] for p in polys]
        psd(gram)
        self.assertTrue(all(v == 0 for v in gram[0]))
        # A physical trace cannot remove an erroneously populated forced-null
        # row by a mixture of weight strictly below one.
        for t in [F(0), F(1, 2), F(999999, 1000000)]:
            with self.assertRaises(ValueError):
                psd([[F(0), (1-t)/100], [(1-t)/100, F(1)]])

    def test_invalid_domain_refused(self):
        with self.assertRaises(ValueError):
            singlet_trace(mono(()), 8, 3)
        with self.assertRaises(ValueError):
            singlet_trace(mono(((1, 8), (0, 8))), 8, 4)

    def test_against_independent_small_bit_state_spin_projector(self):
        from research.transfer_solver_20260915.independent_oracle import matrix
        # Explicit enumeration is confined to this 15-state validation.
        _, s2 = matrix({'modes': 6, 'particles': 2, 'hamiltonian': encode(spin_squared(6))})
        p0 = [[F(i == j)-s2[i][j]/2 for j in range(15)] for i in range(15)]
        operators = [mono(((1, 0), (0, 0))), mono(((1, 0), (1, 3), (0, 0), (0, 3))),
            mono(((1, 0), (1, 3), (0, 1), (0, 2)))]
        for p in operators:
            p = add(p, canonical(adj(p)))
            _, a = matrix({'modes': 6, 'particles': 2, 'hamiltonian': encode(p)})
            exact = sum(p0[i][j]*a[j][i] for i in range(15) for j in range(15))/6
            self.assertEqual(singlet_trace(p, 6, 2), exact)


if __name__ == '__main__':
    unittest.main()
