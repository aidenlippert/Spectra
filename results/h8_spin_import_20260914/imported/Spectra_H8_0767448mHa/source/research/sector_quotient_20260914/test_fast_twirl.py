from fractions import Fraction as F
from itertools import combinations
import unittest
from research.sector_quotient_20260914.fast_twirl import twirl
from research.certificate_scaling.spin_twirl import twirl as original


class FastProjectionTests(unittest.TestCase):
    def test_all_662_number_conserving_six_mode_matrix_units(self):
        count = 0
        for k in range(4):
            for I in combinations(range(6), k):
                for J in combinations(range(6), k):
                    w = tuple((1, i) for i in I)+tuple((0, j) for j in J)
                    p = {w: F(7, 11)}
                    self.assertEqual(twirl(p), original(p))
                    count += 1
        self.assertEqual(count, 662)

    def test_sixteen_mode_sparse_mixture_and_refusal(self):
        p = {((1, 0), (1, 4), (1, 13), (0, 2), (0, 8), (0, 15)): F(2, 9),
             ((1, 0), (1, 15), (0, 3), (0, 14)): F(-3, 7), (): F(5)}
        self.assertEqual(twirl(p), original(p))
        with self.assertRaises(ValueError): twirl({((0, 0), (1, 1)): F(1)})
        with self.assertRaises(ValueError): twirl({(): .5})


if __name__ == '__main__': unittest.main()
