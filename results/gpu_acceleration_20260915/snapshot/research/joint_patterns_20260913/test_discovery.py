"""Numerical-map checks against exact CAR columns and rounded proof replay."""
from fractions import Fraction as F
import unittest
import numpy as np

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import add, canonical, mono
from research.certificate_scaling.polynomial_gram_contraction import contract
from research.joint_patterns_20260913.core import anticommutator, prepare_anticommutators, replay
from research.joint_patterns_20260913.discovery import export_compact
from research.joint_patterns_20260913 import test_core as fixtures


class ProposalMaps(unittest.TestCase):
    def test_contracted_cross_columns_with_and_without_basis_change(self):
        polys = [add(mono(((0, 0),)), mono(((1, 1), (0, 0), (0, 2)), F(2, 3))),
            add(mono(((0, 3),)), mono(((1, 2), (0, 1), (0, 3)), F(-3, 5)))]
        prepared, words, _ = prepare_anticommutators([{'polynomials': polys}])
        rows = sorted(words); lookup = {w: i for i, w in enumerate(rows)}
        for W in (None, np.array([[1., 2.], [-3., 1.]])):
            matrix, ix, _ = contract(prepared[0], polys, lookup, transform=W)
            for z in (np.array([1., -2.]), np.array([3., 4.])):
                coefficients = z if W is None else z@W
                q = add(*({w: F(float(c))*a for w, a in p.items()} for p, c in zip(polys, coefficients)))
                exact = anticommutator(q, q)
                self.assertTrue(np.allclose(matrix@(np.outer(z, z).ravel()[ix]),
                    [float(exact.get(w, 0)) for w in rows], atol=1e-12, rtol=1e-12))

    def test_export_ties_survive_eigenvalue_clipping_and_rounding(self):
        fixtures.JointPatterns.setUpClass(); case = fixtures.JointPatterns
        refs = [[-1, 0], [0, 1], [1, 0]]
        groups = [{'name': 'test', 'kind': 'anticommutator', 'generators': refs},
            {'name': 'base', 'kind': 'square', 'polynomials': [mono(((0, 0),))]}]
        grams = [np.array([[1., .1, .2], [.1, 1., 0.], [.2, 0., -.001]]), np.array([[.25]])]
        cert, clipped = export_compact(case.data, case.tail, groups, grams, [np.eye(3), None], [0.], [])
        self.assertGreater(clipped, 0)
        receipt = replay(case.data, case.tail, cert)
        self.assertLessEqual(receipt['retained']['residual_max_degree'], 4)


if __name__ == '__main__':
    unittest.main()
