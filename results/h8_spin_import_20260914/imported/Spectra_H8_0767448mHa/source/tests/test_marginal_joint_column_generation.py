import contextlib
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_joint_coefficient_constructor import construct as prepare
from experiments.marginal_joint_spinflip import reduce_cached
from experiments.marginal_joint_column_generation import construct


class ColumnGenerationTests(unittest.TestCase):
    def test_small_exact_export_and_input_refusals(self):
        c = build(2, 4, F(1, 3)); data = {k: c[k] for k in ('modes', 'particles', 'hamiltonian')}
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()), \
             patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')), \
             patch('experiments.marginal_spin_constructor.spin_states', side_effect=AssertionError('No states')), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals', side_effect=AssertionError('No lift')):
            root = Path(folder); source = root/'h.json'; source.write_text(json.dumps(data))
            prepare(data, 0, root/'cache', time_limit=10)
            reduce_cached(root/'cache', source, root/'reduced')
            result = construct(data, root/'reduced', root/'proof', time_limit=10)
            self.assertTrue(result['exact_accepted'])
            self.assertTrue((root/'proof/active_columns.json').exists())
            with self.assertRaises(ValueError):
                construct(dict(data, hamiltonian=[]), root/'reduced', root/'bad')
            with self.assertRaises(ValueError):
                construct(data, root/'reduced', root/'bad', batch=0)


if __name__ == '__main__':
    unittest.main()
