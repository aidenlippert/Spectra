import unittest
import numpy as np
from research.acceptance_channels_20260915.accepted_fit import l1_fit


class AcceptedFitTest(unittest.TestCase):
    def test_joint_scalar_and_positive_square_fit_has_known_optimum(self):
        # H=diag(2,3), Q=diag(0,q), and scalar b. The exact optimum is b=2,q=1.
        step, receipt = l1_fit(np.array([[1., 0.], [1., 1.]]),
            np.array([2., 3.]), np.ones(2), np.array([1., 0.]), [None, 0.], 5)
        self.assertIsNotNone(step)
        b, q = step
        lower = b - np.abs(np.array([2., 3.])-np.array([b, b+q])).sum()
        self.assertAlmostEqual(lower, 2., places=8)
        self.assertGreaterEqual(q, 0)
        self.assertFalse(receipt['exact_family_obstruction'])

    def test_weighted_l1_fit_differs_from_least_squares(self):
        step, _ = l1_fit(np.ones((3, 1)), np.array([0., 0., 10.]),
            np.ones(3), np.zeros(1), [None], 5)
        self.assertAlmostEqual(step[0], 0., places=8)

    def test_fixed_rows_remain_in_objective(self):
        _, receipt = l1_fit(np.array([[1.], [0.]]), np.array([0., .25]),
            np.array([1., 2.]), np.zeros(1), [None], 5)
        self.assertAlmostEqual(receipt['objective_Ha'], .5)

    def test_invalid_weights_refused(self):
        with self.assertRaises(ValueError):
            l1_fit(np.ones((1, 1)), np.ones(1), np.array([-1.]), np.zeros(1), [None], 5)


if __name__ == '__main__':
    unittest.main()
