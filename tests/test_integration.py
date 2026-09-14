import json
import unittest
from fractions import Fraction
from itertools import product
import numpy as np
from experiments.pauli import multiply, commutator_i, generator_matrix
from experiments.run import dense_pauli, dense_hamiltonian, integrate


class IntegrationTests(unittest.TestCase):
    def test_exact_algebra_against_independent_dense_matrices(self):
        labels = [''.join(x) for x in product('IXYZ', repeat=2)]
        for p in labels:
            for q in labels:
                phase, r = multiply(p, q)
                np.testing.assert_allclose(dense_pauli(p) @ dense_pauli(q), phase * dense_pauli(r), atol=0)
                coeffs = commutator_i({p: Fraction(1, 3)}, q)
                expected = 1j / 3 * (dense_pauli(p) @ dense_pauli(q) - dense_pauli(q) @ dense_pauli(p))
                actual = dense_hamiltonian(coeffs) if coeffs else np.zeros((4, 4))
                np.testing.assert_allclose(actual, expected, atol=1e-15)

    def test_withheld_coupling_has_explicit_omission(self):
        _, _, res = generator_matrix({'ZXI': 1}, ['ZII', 'YII'])
        self.assertEqual(res['YII'], 2)
        self.assertEqual(commutator_i({'ZXI': 1}, 'YII'), {'XXI': Fraction(2)})

    def test_end_to_end_certificates_and_dynamics(self):
        result = integrate(7)
        json.dumps(result)
        for stage in result['stages']:
            for context in stage['contexts']:
                self.assertTrue(context['energy']['independent_rational_check'])
                self.assertTrue(context['energy']['oracle_numeric_bracket_check'])
                self.assertTrue(context['representation']['exactly_closed_for_each_control'])
                self.assertTrue(context['representation']['trajectory']['numerically_within_bound'])
                self.assertTrue(context['representation']['budget_two_abstained'])
                self.assertTrue(context['oracle_trajectory']['numerically_within_bound'])
                if stage['stage']:
                    self.assertTrue(context['representation']['prior_basis_failure_witnesses'])


if __name__ == '__main__':
    unittest.main()
