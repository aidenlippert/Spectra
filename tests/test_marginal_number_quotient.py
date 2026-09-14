from fractions import Fraction as F
from itertools import combinations
from math import comb
from types import SimpleNamespace
import random
import unittest
from unittest.mock import patch

from experiments.marginal_number_quotient import NumberSliceQuotient, complete_bounded_number_ideals


def add(a, b):
    c = dict(a)
    for k, v in b.items():
        c[k] = c.get(k, F(0))+v
    return {k: v for k, v in c.items() if v}


def multiply(a, b):
    c = {}
    for i, x in a.items():
        for j, y in b.items():
            c[i | j] = c.get(i | j, F(0))+x*y
    return {k: v for k, v in c.items() if v}


def evaluate(poly, assignment):
    return sum((v for k, v in poly.items() if k & assignment == k), F(0))


class NumberQuotientTests(unittest.TestCase):
    def test_exact_witness_and_degree_on_every_small_slice(self):
        rng = random.Random(1729)
        for m, p, d in [(4, 2, 2), (6, 3, 3), (8, 4, 2), (6, 2, 2), (6, 4, 2)]:
            q = NumberSliceQuotient(list(range(m)), p, d)
            self.assertEqual(q.stats['quotient_dimension'], comb(m, d))
            polynomial = {s: F(rng.randrange(-5, 6), 7) for s in q.monomials}
            polynomial = {k: v for k, v in polynomial.items() if v}
            reduced, ideal = q.reduce(polynomial)
            shift = {0: -p, **{1 << i: 1 for i in range(m)}}
            self.assertEqual(add(reduced, multiply(shift, ideal)), polynomial)
            self.assertLessEqual(max(map(int.bit_count, reduced), default=0), d)
            self.assertLessEqual(max(map(int.bit_count, ideal), default=0), d-1)
            for chosen in combinations(range(m), p):
                state = sum(1 << i for i in chosen)
                self.assertEqual(evaluate(polynomial, state), evaluate(reduced, state))
            for support in q.monomials:
                if support.bit_count() < d:
                    zero = multiply(shift, {support: F(2, 3)})
                    self.assertEqual(q.reduce(zero)[0], {})

    def test_two_spin_exact_identity_without_actions_or_full_population_lift(self):
        ring = SimpleNamespace(modes=16, target=4,
            spin_masks=[sum(1 << i for i in range(s, 16, 2)) for s in (0, 1)])
        polynomial = {0: F(3, 7), 1: F(-2, 3), 3: F(4, 9), 21: F(7, 11),
                      42: F(2, 13), 51: F(-3, 17), 85: F(1, 19)}
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError):
            reduced, ideals, receipts = complete_bounded_number_ideals(ring, polynomial)
        reconstructed = dict(reduced)
        for spin, ideal in enumerate(ideals):
            shift = {0: -4, **{1 << i: 1 for i in range(spin, 16, 2)}}
            reconstructed = add(reconstructed, multiply(shift, ideal))
        self.assertEqual(reconstructed, polynomial)
        self.assertLessEqual(max(map(int.bit_count, reduced), default=0), 4)
        self.assertTrue(all(x['configuration_evaluations'] == 0 for x in receipts))
        self.assertLess(receipts[1]['maximum_degree'], ring.target)

    def test_boundaries_and_refusals(self):
        q = NumberSliceQuotient([0, 2, 4, 6], 2, 0)
        self.assertEqual(q.reduce({0: F(3, 5)}), ({0: F(3, 5)}, {}))
        for args in [([], 1, 1), ([0, 0], 1, 1), ([True, 1], 1, 1),
                     ([0, 1], 0, 0), ([0, 1, 2, 3], 2, 3),
                     ([0, 1, 2, 3], 2, 2, 2)]:
            with self.assertRaises(ValueError):
                NumberSliceQuotient(*args)
        for polynomial in [{1: 1}, {-1: 1}, {True: 1}]:
            with self.assertRaises(ValueError):
                q.reduce(polynomial)


if __name__ == '__main__':
    unittest.main()
