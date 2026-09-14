from fractions import Fraction as F
from math import prod
from itertools import combinations, product
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from experiments.marginal_joint_coefficient_constructor import feature_orbits, vanishing_feature, sector_mean, prepare, construct, solve_prepared, export
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_polynomial_metric import JointPolynomial


class JointCoefficientConstructorTests(unittest.TestCase):
    def test_generated_orbits_and_vanishing_feature_identity(self):
        orbits = feature_orbits(4)
        self.assertEqual(len({p for orbit in orbits for p in orbit}), sum(1 for p in product(range(3), repeat=4) if sum(p) <= 4))
        for orbit in orbits:
            self.assertEqual(set(orbit), set(p[::-1] for p in orbit))
            metric = vanishing_feature(orbit, 4)
            self.assertTrue(all(sum(t['powers']) <= 6 and any(t['powers']) for t in metric['terms']))
            for q in [(-1, 0, 0, 1), (1, 1, -1, -1), (0, 0, 0, 0)]:
                value = lambda p: prod(x**k for x, k in zip(q, p))
                actual = sum(F(t['coefficient'], 2)*value(t['powers']) for t in metric['terms'])
                expected = F(sum(x*x for x in q), 2)*sum(value(p) for p in orbit)
                self.assertEqual(actual, expected)

    def test_analytic_mean_matches_small_sector_and_Q_normalization(self):
        sites, p = 4, 2
        states = [sum(1 << (2*i) for i in a)+sum(1 << (2*i+1) for i in b)
                  for a in combinations(range(sites), p) for b in combinations(range(sites), p)]
        for mask in [0, 1, 3, 5, 21, 85, 15]:
            exact = F(sum(mask&s == mask for s in states), len(states))
            self.assertEqual(sector_mean({mask: F(1)}, sites, p), exact)
        c = build(4, 4, F(1, 3)); ring = JointPolynomial(c)
        q = [{0: -1, 1 << (2*i): 1, 1 << (2*i+1): 1} for i in range(sites)]
        w = {m: F(v, 2) for m, v in ring.metric_polynomial(q).items()}
        self.assertEqual(F(6, 5)*sector_mean(w, sites, p), F(6, 5))

    def test_small_bare_H_end_to_end_and_input_refusals(self):
        c = build(2, 4, F(1, 3)); data = {k: c[k] for k in ['modes', 'particles', 'hamiltonian']}
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()), \
             patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
             patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state list')), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No lifting')):
            out = Path(folder)
            result = construct(data, F(0), out, time_limit=10)
            self.assertTrue(result['success'])
            with patch('experiments.marginal_joint_coefficient_constructor.prepare', side_effect=AssertionError('Reuse prepared coefficients')):
                resumed = solve_prepared(data, F(0), out, time_limit=10, method='highs-ds')
                self.assertTrue(resumed['success'])
            with self.assertRaises(ValueError):
                solve_prepared(data, F(1), out, time_limit=10)
            proposal = json.loads((out/'proposal.json').read_text())
            proof = export(data, proposal, out/'proof')
            with self.assertRaises(ValueError):
                export(dict(data, hamiltonian=[]), proposal, out/'invalid')
            self.assertGreater(F(proof['weight_positivity']['lower']), 0)
            self.assertGreater(F(proof['numerator_positivity']['lower']), 0)
            self.assertEqual(proof['compilation']['determinant_actions'], 0)
            with self.assertRaises(ValueError):
                prepare(dict(data, polynomial_metric=c['polynomial_metric']), F(0))
            with self.assertRaises(ValueError):
                prepare(dict(data, modes=100), F(0))
            with self.assertRaises(ValueError):
                prepare(data, F(0), max_columns=1)


if __name__ == '__main__':
    unittest.main()
