from fractions import Fraction as F
from itertools import product
import unittest
import numpy as np
from experiments.v4_moments import signed_gram_plan
from experiments.v4_planner import action_values_h2


class MomentTests(unittest.TestCase):
    def test_exact_perfect_and_uninformative(self):
        p, s = (F(1, 2), F(1, 2)), (-1, 1)
        table = ((F(0), F(1, 2)), (F(1), F(1, 2)))
        result = signed_gram_plan(p, table, s)
        self.assertEqual(result['values'], (F(0), F(0)))
        self.assertEqual(result['work']['G_products'], 6)
        empty = ((F(1, 2),), (F(1, 2),))
        self.assertEqual(signed_gram_plan(p, empty, s)['values'], (F(1, 2),))

    def test_first_signed_moments_are_insufficient(self):
        p, s = (F(1, 4),) * 4, (1, 1, -1, -1)
        world_x = tuple((F(a), F(a)) for a in (1, 0, 1, 0))
        world_y = tuple((F(a), F(b)) for a, b in zip((1, 0, 1, 0), (1, 0, 0, 1)))
        x, y = signed_gram_plan(p, world_x, s), signed_gram_plan(p, world_y, s)
        self.assertEqual(x['moments']['first'], y['moments']['first'])
        self.assertEqual(min(x['values']), F(1, 2))
        self.assertEqual(min(y['values']), F(0))

    def test_random_biased_priors_match_independent_recursion(self):
        rng = np.random.default_rng(47)
        for _ in range(24):
            raw = [int(x) for x in rng.integers(1, 9, size=6)]
            prior = tuple(F(x, sum(raw)) for x in raw)
            table = tuple(tuple(F(int(x), 6) for x in row) for row in rng.integers(0, 7, size=(6, 5)))
            signs = (-1, 1, -1, 1, -1, 1)
            result = signed_gram_plan(prior, table, signs)
            self.assertEqual(result['values'], action_values_h2(prior, table, signs))

    def test_invalid_exact_inputs(self):
        with self.assertRaises(ValueError): signed_gram_plan((1.,), ((F(1, 2),),), (1,))
        with self.assertRaises(ValueError): signed_gram_plan((F(1),), ((),), (1,))


if __name__ == '__main__': unittest.main()
