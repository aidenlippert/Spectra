"""Individually identical experiment statistics can require opposite choices."""
from fractions import Fraction as F
from itertools import product, combinations_with_replacement
from math import prod
from types import SimpleNamespace
import unittest
import numpy as np
from experiments.v4_run import features
from experiments.v4_moments import signed_gram_plan


def indistinguishable_contexts():
    bits = list(product((0, 1), repeat=3))
    signs = tuple(1 if u == v else -1 for u, v, w in bits)
    prior = (F(1, 8),) * 8
    first = tuple(tuple(map(F, (u, v, w, w))) for u, v, w in bits)
    second = tuple(tuple(map(F, (w, w, u, v))) for u, v, w in bits)
    common = {'prior': prior, 'signs': signs, 'n': 3,
              'actions': [SimpleNamespace(prep='ZII')] * 4}
    return dict(common, table=first), dict(common, table=second)


class FeatureObstructionTests(unittest.TestCase):
    def test_all_twelve_features_equal_but_optima_disjoint(self):
        left, right = indistinguishable_contexts()
        x, _ = features(left); y, _ = features(right)
        np.testing.assert_array_equal(x, y)
        a = signed_gram_plan(left['prior'], left['table'], left['signs'])
        b = signed_gram_plan(right['prior'], right['table'], right['signs'])
        self.assertEqual(a['values'], (F(0), F(0), F(1, 2), F(1, 2)))
        self.assertEqual(b['values'], (F(1, 2), F(1, 2), F(0), F(0)))
        # Every fixed action, and therefore every randomized action independent
        # of the hidden world, has exactly 1/4 mean regret across these worlds.
        self.assertTrue(all((v + w) / 2 == F(1, 4) for v, w in zip(a['values'], b['values'])))


    def test_higher_moment_hierarchy_by_enumeration(self):
        for h in range(2, 6):
            worlds = list(product((0, 1), repeat=h + 1))
            labels = ([(-1) ** sum(row[:h]) for row in worlds],
                      [(-1) ** row[h] for row in worlds])
            for order in range(h):
                for indices in combinations_with_replacement(range(h), order):
                    for signs in labels:
                        moment = sum(s * prod(row[i] for i in indices)
                                     for s, row in zip(signs, worlds))
                        self.assertEqual(moment, 0)
            risks = []
            for signs in labels:
                groups = {}
                for row, sign in zip(worlds, signs):
                    masses = groups.setdefault(row[:h], [0, 0])
                    masses[int(sign == 1)] += 1
                risks.append(F(sum(min(masses) for masses in groups.values()), len(worlds)))
            self.assertEqual(risks, [F(0), F(1, 2)])


if __name__ == '__main__': unittest.main()
