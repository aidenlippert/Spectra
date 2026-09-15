import unittest
from fractions import Fraction
from experiments.v7_identifiability import identify_query, verify_row_certificate, verify_nullspace_witness, verify_compatible_models


class IdentifiabilityTests(unittest.TestCase):
    def test_exact_identified_query_rank_deficient(self):
        X = [[1, 0, 1], [0, 1, 1]]
        c = identify_query(X, [1, 1, 2], y=[3, 4])
        self.assertTrue(c.identified)
        self.assertTrue(verify_row_certificate(X, [1, 1, 2], c.weights))

    def test_unidentified_query_has_compatible_models(self):
        c = identify_query([[1, 0], [2, 0]], [0, 1], y=[3, 6], margin=2)
        self.assertFalse(c.identified)
        self.assertTrue(verify_nullspace_witness([[1, 0], [2, 0]], c.nullspace_witness, [0, 1]))
        self.assertEqual(sum(c.theta_plus[i]*[0, 1][i] for i in range(2)) - sum(c.theta_minus[i]*[0, 1][i] for i in range(2)), Fraction(2))

    def test_corrupt_witness_is_rejected(self):
        self.assertFalse(verify_nullspace_witness([[1, 0]], [1, 1], [0, 1]))
        self.assertFalse(verify_row_certificate([[1, 0]], [1, 1], [1]))

    def test_noise_compatibility_and_identified_bound(self):
        X = [[1, 0], [0, 1]]
        c = identify_query(X, [1, 2], y=[3, 4], theta0=[3, 4], eta=Fraction(1, 10))
        self.assertTrue(c.identified)
        self.assertEqual(c.error_bound, Fraction(3, 10))
        with self.assertRaises(ValueError): identify_query([[1, 0]], [0, 1], y=[3], theta0=[9, 0], eta=1)

    def test_incompatible_exact_data_refused(self):
        with self.assertRaises(ValueError): identify_query([[1, 0], [1, 0]], [0, 1], y=[1, 2])

    def test_identified_incompatible_and_bad_arguments(self):
        with self.assertRaises(ValueError): identify_query([[1, 0], [1, 0]], [1, 0], y=[1, 2])
        with self.assertRaises(ValueError): identify_query([[1, 0]], [0, 1], eta=-1, theta0=[0, 0], y=[0])
        with self.assertRaises(ValueError): identify_query([[1, 0]], [0, 1], margin=0)
        self.assertFalse(verify_row_certificate([[1, 0]], [1, 2], [1]))
        self.assertFalse(verify_nullspace_witness([[1, 0]], [0], [0, 1]))
        self.assertFalse(verify_nullspace_witness([[1, 0]], [0, 1], [0]))
        self.assertFalse(verify_nullspace_witness([[1, 0]], [0, 1], [0, 1, 2]))
        self.assertFalse(verify_row_certificate([[1, 0]], ["bad", 1], [1]))
        self.assertFalse(verify_nullspace_witness([[1, 0]], [0, 1], ["bad", 1]))
        self.assertFalse(verify_compatible_models([[1, 0]], [0], [0, 1], [0], [0], 1, 2))

    def test_unchanged_input_identified_query(self):
        self.assertTrue(identify_query([[1, 50, 1], [2, 50, 1]], [3, 50, 1]).identified)
        self.assertFalse(identify_query([[1, 50, 1], [2, 50, 1]], [3, 75, 1]).identified)


if __name__ == '__main__': unittest.main()
