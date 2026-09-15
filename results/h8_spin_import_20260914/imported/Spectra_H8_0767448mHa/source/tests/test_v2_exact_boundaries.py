from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import unittest
from experiments.v2_laws import CalibrationRecord, learn, verify_certificate, required_shots, _exact_union_bound
from experiments.v2_policy import certificate, verify


class ExactBoundaryTests(unittest.TestCase):
    def test_binomial_tail_against_exhaustive_noise_sequences(self):
        for k in (1, 3, 5, 7):
            for eta in (F(0), F(1, 10), F(2, 5)):
                total = F(0)
                for bits in product((0, 1), repeat=k):
                    if sum(bits) > k // 2:
                        total += eta ** sum(bits) * (1 - eta) ** (k - sum(bits))
                self.assertEqual(_exact_union_bound(3, k, eta), 3 * total)

    def test_all_failed_map_complementarity_cases(self):
        v = F(1, 2)
        for ca, cb in product((0, 1), repeat=2):
            r0 = F(1, 2) - v / 8
            ra = F(1, 2) - v * ca / 4
            rb = F(1, 2) - v * cb / 4
            rab = F(1, 2) - v * ca * cb / 2
            j = ra + rb - r0 - rab
            self.assertEqual(j, v / 8 if ca == cb else -v / 8)
            self.assertGreaterEqual(ra, F(3, 8))
            self.assertGreaterEqual(rb, F(3, 8))
            self.assertLessEqual(rab, F(1, 2))

    def test_malformed_learning_certificates_and_horizons(self):
        records = (CalibrationRecord((1,), (1,)),)
        good = learn(records, F(0), F(1, 100))
        for change in ({'mask': (1.0,)}, {'mask': (True,)}, {'xor_operations': 99},
                       {'exact_error_bound': '1/7'}, {'reason': 'unsupported'}, {'shots_per_context': 3}):
            self.assertFalse(verify_certificate(records, replace(good, **change), F(0), F(1, 100)))
        for h in (1.5, True, '1', -1, 6):
            with self.assertRaises(ValueError): certificate('none', h)
        for d in (True, 0, 257):
            with self.assertRaises(ValueError): required_shots(d, F(1, 10), F(1, 100))


if __name__ == '__main__': unittest.main()
