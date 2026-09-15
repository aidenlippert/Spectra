import unittest
from fractions import Fraction as F
from experiments.v7_headroom import model, TOL
from experiments.v7_certificate import Generator, check_certificate, BudgetExceeded
from experiments.v7_adaptive_taylor import adaptive_taylor
from experiments.v8_integer_taylor import fraction_free_taylor


class FractionFreeTaylorTests(unittest.TestCase):
    def test_matches_adaptive_and_checks(self):
        for n, gamma, duration in ((3, F(0), F(1, 5)), (4, F(2), F(1, 2))):
            h, initial = model(n, 'xxz')
            old_piece, old_witness, _ = adaptive_taylor(Generator(h, gamma, n, max_terms=512), initial, duration, TOL)
            piece, witness, _ = fraction_free_taylor(Generator(h, gamma, n, max_terms=512), initial, duration, TOL)
            self.assertEqual(piece, old_piece)
            self.assertEqual(witness['claimed_bound'], old_witness['claimed_bound'])
            result = check_certificate(Generator(h, gamma, n, max_terms=512), initial, [piece], witness, TOL, expected_time=duration)
            self.assertEqual(result['status'], 'certified')

    def test_stationary_rational_initial(self):
        piece, witness, _ = fraction_free_taylor(Generator({}, F(0), 1), {'X': F(1, 3)}, F(1), F(0))
        self.assertEqual(piece.coefficients, ({'X': F(1, 3)},))
        self.assertEqual(witness['claimed_bound'], '0')

    def test_cache_and_live_refusals(self):
        h, initial = model(6, 'xxz')
        with self.assertRaises(BudgetExceeded):
            fraction_free_taylor(Generator(h, F(0), 6, max_terms=512), initial, F(1, 2), TOL)
        with self.assertRaises(BudgetExceeded):
            fraction_free_taylor(Generator(h, F(2), 6, max_terms=512), initial, F(1, 5), TOL)

    def test_invalid_labels_rejected(self):
        with self.assertRaises(ValueError):
            fraction_free_taylor(Generator({}, F(0), 1), {'Q': F(1)}, F(1), F(0))


if __name__ == '__main__':
    unittest.main()
