import contextlib
from fractions import Fraction as F
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_fixed_metric_pricing import solve_part, combine, IncrementalMaster


class FixedMetricPricingTests(unittest.TestCase):
    def test_native_append_preserves_feasibility_and_removes_residual(self):
        import numpy as np
        from scipy.sparse import csc_matrix
        master=IncrementalMaster(np.array([1.,1.]))
        first=master.run(csc_matrix([[1.],[0.]]),10)
        self.assertTrue(first.success);self.assertAlmostEqual(first.fun,1.)
        second=master.run(csc_matrix([[0.],[1.]]),10)
        self.assertTrue(second.success);self.assertAlmostEqual(second.fun,0.)
        self.assertTrue(np.allclose(second.x[:2],[1.,1.]))
        self.assertEqual(master.columns,2)

    def test_independent_parts_full_replay_and_binding(self):
        raw=build(2,4,F(1,3));candidate={k:raw[k] for k in ('modes','particles','hamiltonian','target_lower','polynomial_metric')}
        with tempfile.TemporaryDirectory() as folder, contextlib.redirect_stdout(io.StringIO()), \
             patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No states')), \
             patch('experiments.marginal_spin_constructor.spin_states',side_effect=AssertionError('No states')), \
             patch('experiments.marginal_polynomial_metric.complete_number_ideals',side_effect=AssertionError('No full lift')):
            out=Path(folder)
            for name in ('weight','numerator'):
                self.assertTrue(solve_part(candidate,name,out/name,time_limit=10,native=True)['exact_accepted'])
            receipt=combine(candidate,out)
            self.assertEqual(receipt['complement_lower'],candidate['target_lower'])
            self.assertGreater(F(receipt['weight_positivity']['lower']),0)
            self.assertGreater(F(receipt['numerator_positivity']['lower']),0)
            with self.assertRaises(ValueError):combine(dict(candidate,target_lower='0'),out)
            with self.assertRaises(ValueError):solve_part(candidate,'unknown',out)
            with self.assertRaises(ValueError):solve_part(candidate,'weight',out,batch=0)
            with self.assertRaises(ValueError):solve_part(dict(candidate,inherited_proof={}), 'weight', out)


if __name__=='__main__':unittest.main()
