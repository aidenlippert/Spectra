from fractions import Fraction as F
from itertools import combinations
import random
import unittest
from research.interacting_scaling_20260915.spin_patterns import pattern, twirl
from research.sector_quotient_20260914.fast_twirl import twirl as reference


class SpinPatternTests(unittest.TestCase):
    def test_all_small_balanced_words_and_large_relabelings(self):
        pattern.cache_clear()
        words = []
        for degree in range(4):
            for left in combinations(range(6), degree):
                for right in combinations(range(6), degree):
                    words.append(tuple((1, i) for i in left)+tuple((0, i) for i in right))
        rng = random.Random(325)
        for _ in range(150):
            left = sorted(rng.sample(range(24), 3)); right = sorted(rng.sample(range(24), 3))
            words.append(tuple((1, i) for i in left)+tuple((0, i) for i in right))
        for w in words:
            self.assertEqual(twirl({w: F(7, 13)}), reference({w: F(7, 13)}))
        self.assertGreater(pattern.cache_info().hits, 0)
        self.assertEqual(twirl({words[0]: F(1, 7), words[-1]: F(-3, 2)}),
            reference({words[0]: F(1, 7), words[-1]: F(-3, 2)}))

    def test_invalid_inputs_keep_the_original_refusals(self):
        for p in [{((0, 0),): F(1)}, {((0, 0), (1, 0)): F(1)},
                  {((1, 2), (1, 0)): F(1)}, {(): .5}]:
            with self.assertRaises(ValueError): twirl(p)


if __name__ == '__main__':
    unittest.main()
