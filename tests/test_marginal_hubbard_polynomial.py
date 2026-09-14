import copy
from fractions import Fraction as F
from itertools import combinations
import unittest
from unittest.mock import patch

from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_polynomial_metric import JointPolynomial, replay
from experiments.marginal_spin_reduction import SpinZeroOracle


class AnalyticHubbardCertificateTest(unittest.TestCase):
    def test_exact_size_transfer_without_actions_or_sector_lifting(self):
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                   side_effect=AssertionError('No physical rows permitted')), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals',
                   side_effect=AssertionError('No sector-sized lifting permitted')), \
             patch('experiments.marginal_spin_constructor.spin_states',
                   side_effect=AssertionError('No sector enumeration permitted')):
            for sites in (2, 4, 6, 8, 10):
                c = build(sites, 4, F(1, 3))
                r = replay(c)
                self.assertEqual(r['weight_positivity']['residual_l1'], '0')
                self.assertEqual(r['numerator_positivity']['residual_l1'], '0')
                self.assertEqual(F(r['weight_positivity']['lower']), 1)
                self.assertEqual(F(r['numerator_positivity']['lower']), F(1, 1000))
                self.assertTrue(r['compilation']['projector_free_zero_metric'])
                self.assertLessEqual(r['compilation']['numerator_degree'], 4)
                self.assertLessEqual(len(c['numerator_proof']['positive_indicators']), 4*sites*(sites-1)+12*(sites-1))
                self.assertTrue(all(a['required'].bit_count() <= 4 for a in c['numerator_proof']['positive_indicators']))
                self.assertTrue(all(a['mask'].bit_count() <= 2 for v in c['numerator_proof']['number_multipliers'] for a in v))

    def test_compiled_identity_against_actual_car_rows_including_valence(self):
        c = build(4, 2, F(1, 3))
        oracle = SpinZeroOracle(c)
        ring = JointPolynomial(c)
        _, polynomial, cost = ring.compile(F(c['target_lower']))
        for alpha in combinations(range(4), 2):
            for beta in combinations(range(4), 2):
                s = sum(1 << (2*i) for i in alpha)+sum(1 << (2*i+1) for i in beta)
                d = sum(((s >> (2*i)) & 3) == 3 for i in range(4))
                row = oracle.action(s)
                expected = (row.get(s, F(0))-F(c['target_lower']))*d
                expected -= sum((abs(v)*sum(((t >> (2*i)) & 3) == 3 for i in range(4))
                                 for t, v in row.items() if t != s), F(0))
                actual = sum((F(v, cost['numerator_scale']) for mask, v in polynomial.items() if mask & s == mask), F(0))
                self.assertEqual(actual, expected)
                if d:
                    self.assertGreaterEqual(actual, F(1, 1000))

    def test_boundaries_and_corrupted_claims(self):
        for U, t in [(0, 0), (0, 1), (1, 0)]:
            r = replay(build(4, U, t))
            self.assertEqual(r['numerator_positivity']['residual_l1'], '0')
        for args in [(True, 1, 1), (3, 1, 1), (34, 1, 1), (4, -1, 1),
                     (4, 1, -1), (4, 1., 1), (4, 1, 1, 0),
                     (4, 1, 1, F(1, 10**17))]:
            with self.assertRaises(ValueError):
                build(*args)
        source = build(4, 4, F(1, 3))
        bad = copy.deepcopy(source)
        bad['target_lower'] = str(F(bad['target_lower'])+1)
        with self.assertRaises(ValueError):
            replay(bad)
        bad = copy.deepcopy(source)
        bad['hamiltonian'] = build(4, 3, F(1, 3))['hamiltonian']
        with self.assertRaises(ValueError):
            replay(bad)
        bad = copy.deepcopy(source)
        bad['numerator_proof']['positive_indicators'][0]['weight'] = -1
        with self.assertRaises(ValueError):
            replay(bad)


if __name__ == '__main__':
    unittest.main()
