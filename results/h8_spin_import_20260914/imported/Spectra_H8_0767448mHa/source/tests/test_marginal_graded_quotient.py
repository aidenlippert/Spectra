from fractions import Fraction as F
from itertools import combinations
import contextlib
import io
import unittest
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_polynomial_metric import JointPolynomial
from experiments.marginal_joint_coefficient_constructor import CoefficientQuotient, prepare_features
from experiments.marginal_moment_pricing import MomentDictionary


class GradedQuotientTests(unittest.TestCase):
    def ring(self, sites):
        return JointPolynomial(build(sites, 4, F(1, 3)))

    def test_filtered_quotient_preserves_small_slice_polynomials(self):
        ring = self.ring(4); q = CoefficientQuotient(ring, 2)
        self.assertEqual(len(q.basis), 20)
        states = [sum(1 << (2*i) for i in a)+sum(1 << (2*i+1) for i in b)
                  for a in combinations(range(4), 2) for b in combinations(range(4), 2)]
        for degree in range(3):
            polynomial = {sum(1 << i for i in bits): F(j+1, 7)
                          for j, bits in enumerate(combinations(range(8), degree))}
            reduced = q.normal(polynomial)
            self.assertTrue(all(m in q.index and m.bit_count() <= degree for m in reduced))
            for state in states:
                evaluate = lambda p: sum(v for m, v in p.items() if state&m == m)
                self.assertEqual(evaluate(polynomial), evaluate(reduced))

    def test_eight_site_rank_number_relations_and_degree_refusals(self):
        ring = self.ring(8); q = CoefficientQuotient(ring, 6)
        self.assertEqual(len(q.basis), 3920)
        self.assertEqual(q.total_degree, 6)
        for support in (0, 5, 31, 341):
            for spin in ring.spin_masks:
                polynomial = ring.multiply({support: F(1)}, {0: -ring.target, **{1 << i: 1 for i in range(ring.modes) if spin&(1 << i)}})
                self.assertEqual(q.normal(polynomial), {})
        with self.assertRaisesRegex(ValueError, 'total-degree'):
            q.normal({127: F(1)})
        for degree, budget in [(-1, 4096), (6, 3919), (6, 0)]:
            with self.assertRaises(ValueError):
                CoefficientQuotient(ring, degree, max_rows=budget)

    def test_reduced_metric_family_has_supported_hubbard_numerator(self):
        c = build(8, 4, F(1, 3)); data = {k: c[k] for k in ('modes', 'particles', 'hamiltonian')}
        with contextlib.redirect_stdout(io.StringIO()):
            ring, q, orbits, weight, numerator, means = prepare_features(data, -16, feature_degree=2)
        self.assertEqual(len(q.basis), 3920)
        self.assertTrue(all(m.bit_count() <= 6 for p in weight+numerator for m in p))
        self.assertEqual(len(orbits), len(means))
        with contextlib.redirect_stdout(io.StringIO()):
            dictionary = MomentDictionary(data, -19, feature_degree=0, proof_degree=4)
        self.assertEqual(dictionary.stats['joint_quotient_dimension'], 1260)
        self.assertEqual(dictionary.localizer_degree, 2)
        self.assertTrue(all(label[1].bit_count() <= (2 if label[0] == 'charge' else 4) for label in dictionary.labels()))
        with self.assertRaises(ValueError):
            prepare_features(data, -16, feature_degree=5)
        with self.assertRaises(ValueError):
            prepare_features(data, -16, quotient_degree=7)


if __name__ == '__main__':
    unittest.main()
