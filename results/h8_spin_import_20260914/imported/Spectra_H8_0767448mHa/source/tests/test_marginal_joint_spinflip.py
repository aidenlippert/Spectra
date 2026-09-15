import contextlib
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_joint_spinflip import flip_mask, flip_label, _orbits, reduce_cached, solve
from experiments.marginal_joint_coefficient_constructor import construct
from experiments.marginal_hubbard_polynomial import build


class SpinflipTests(unittest.TestCase):
    def test_mask_and_label_action_with_refusals(self):
        for mask in range(1 << 8):
            self.assertEqual(flip_mask(flip_mask(mask, 4), 4), mask)
        self.assertEqual(flip_label(['positive', 5, 1], 2), ['positive', 10, 2])
        for mask, sites in [(-1, 2), (16, 2), (1, 0), (True, 2)]:
            with self.assertRaises(ValueError):
                flip_mask(mask, sites)
        self.assertEqual(_orbits([1, 0, 2]), [[0, 1], [2]])
        for permutation in [[0, 0], [1, 2, 0], [2, 0]]:
            with self.assertRaises(ValueError):
                _orbits(permutation)

    def test_small_reduced_discovery_full_export_and_cache_refusals(self):
        import numpy as np
        from scipy.sparse import load_npz, save_npz
        c = build(2, 4, F(1, 3))
        data = {k: c[k] for k in ('modes', 'particles', 'hamiltonian')}
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()), \
             patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
             patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No state generation')), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No full-population lift')):
            root = Path(folder); source = root/'hamiltonian.json'; source.write_text(json.dumps(data))
            cached, reduced = root/'cached', root/'reduced'
            self.assertTrue(construct(data, F(0), cached, time_limit=10)['success'])
            meta = reduce_cached(cached, source, reduced)
            self.assertEqual(meta['original_shape'][0], 9)
            self.assertEqual(meta['reduced_shape'][0], 7)
            self.assertEqual(meta['exact_feature_symmetries'], len(meta['orbits']))
            receipt = solve(data, reduced, time_limit=10)
            self.assertTrue(receipt['exact_accepted'])
            self.assertGreater(F(receipt['exact_receipt']['weight_positivity']['lower']), 0)
            self.assertGreater(F(receipt['exact_receipt']['numerator_positivity']['lower']), 0)
            self.assertEqual(receipt['exact_receipt']['compilation']['determinant_actions'], 0)
            proposal = json.loads((reduced/'proposal.json').read_text())
            for name in ('weight', 'numerator'):
                weights = {tuple(label): value for label, value in zip(proposal[name+'_labels'], proposal[name+'_values'])}
                for label, value in weights.items():
                    self.assertEqual(value, weights[tuple(flip_label(label, 2))])
            source.write_text(json.dumps(dict(data, hamiltonian=[])))
            with self.assertRaisesRegex(ValueError, 'bound'):
                reduce_cached(cached, source, reduced)
            source.write_text(json.dumps(data))
            matrix = load_npz(cached/'model.npz').tolil()
            matrix[0, 0] += 1
            save_npz(cached/'model.npz', matrix.tocsc())
            with self.assertRaisesRegex(ValueError, 'exact compilation'):
                reduce_cached(cached, source, reduced)


if __name__ == '__main__':
    unittest.main()
