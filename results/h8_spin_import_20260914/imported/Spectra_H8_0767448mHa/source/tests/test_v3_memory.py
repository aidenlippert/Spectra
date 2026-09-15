from copy import deepcopy
from fractions import Fraction as F
import unittest
import numpy as np
from experiments.v3_memory import certificate, verify


class MemoryTests(unittest.TestCase):
    def test_full_declared_block_and_independent_dense_trajectory(self):
        c = certificate([[0]], [2], [[F(1, 4)]], [1], times=(F(1, 10), 1, 10), retained_rates=[1])
        self.assertTrue(verify(c))
        L = np.array([[-1, .25], [-.25, -2]])
        vals, vec = np.linalg.eig(L)
        for row in c['times']:
            t = float(F(row['time']))
            actual = (vec @ (np.exp(vals * t) * np.linalg.solve(vec, [1, 0])))[0]
            error = abs(actual - np.exp(-t))
            self.assertLessEqual(error, min(float(F(row['finite_error_bound'])), float(F(row['uniform_error_bound']))) + 1e-14)

    def test_tampering_and_missing_damping_rejected(self):
        c = certificate([[0]], [1], [[2]], [1])
        for field, value in [('B_norm_upper', '1'), ('retained_rates', ['0']),
                             ('rates', ['-1']), ('discarded_initial_is_zero', False)]:
            bad = deepcopy(c); bad[field] = value
            self.assertFalse(verify(bad))
        bad = deepcopy(c); bad['times'][0]['time'] = '-1'
        self.assertFalse(verify(bad))
        with self.assertRaises(ValueError): certificate([[0]], [1], [[1]], [1], retained_rates=[0])

    def test_bad_shape_and_skew_rejected(self):
        with self.assertRaises(ValueError): certificate([[1]], [1], [[1]], [1])
        with self.assertRaises(ValueError): certificate([[0]], [1], [[1, 2]], [1])


if __name__ == '__main__': unittest.main()
