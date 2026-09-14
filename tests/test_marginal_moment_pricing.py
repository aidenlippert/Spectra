import contextlib
from fractions import Fraction as F
import io
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_moment_pricing import MomentDictionary, construct


class MomentPricingTests(unittest.TestCase):
    def data(self, sites):
        c = build(sites, 4, F(1, 3))
        return {k: c[k] for k in ('modes', 'particles', 'hamiltonian')}

    def test_pricing_matches_all_direct_atom_columns(self):
        with contextlib.redirect_stdout(io.StringIO()):
            d = MomentDictionary(self.data(4), F(0))
        rng = np.random.default_rng(173)
        dual = rng.normal(size=2*d.qrows+1)
        active = set()
        selected, maximum, checked = d.price(dual, active, 100000)
        expected = set(); expected_max = 0.
        for label in d.labels():
            column = d.column(label)
            for block in (0, 1):
                score = -float(column @ dual[block*d.qrows:(block+1)*d.qrows])
                expected_max = max(expected_max, score)
                if score > 1e-8:
                    expected.add((block, label))
        self.assertEqual(set(selected), expected)
        self.assertAlmostEqual(maximum, expected_max, places=10)
        self.assertEqual(checked, 2*len(d.labels()))
        selected_again, _, _ = d.price(dual, set(selected), 100000)
        self.assertEqual(selected_again, [])
        self.assertFalse(d.stats['materialized_all_atom_matrix'])

    def test_small_bare_online_export_without_full_dictionary(self):
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()), \
             patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
             patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')), \
             patch('experiments.marginal_joint_coefficient_constructor.prepare', side_effect=AssertionError('No full dictionary')):
            result = construct(self.data(2), F(0), folder, time_limit=10)
            self.assertTrue(result['exact_accepted'])
            self.assertFalse(result['materialized_all_atom_matrix'])
            with self.assertRaises(ValueError):
                MomentDictionary(dict(self.data(2), inherited_metric={}), F(0))


if __name__ == '__main__':
    unittest.main()
