import unittest
from fractions import Fraction as F

from experiments.v8_time_reduction_probe import attempt


class TimeReductionWholeHorizonTests(unittest.TestCase):
    def test_segmented_attempt_is_bound_to_full_horizon(self):
        row = attempt(3, F(0), F(1, 5), 2, 8)
        self.assertIn(row['status'], ('certified', 'over_tolerance'))
        self.assertEqual(row['segments'], 2)
        self.assertEqual(row['order'], 8)
        if row['status'] == 'certified':
            self.assertEqual(row['best']['status'], 'certified')
            self.assertLessEqual(F(row['best']['bound']), F(1, 1000))


if __name__ == '__main__':
    unittest.main()
