import unittest
import math
from fractions import Fraction as F
from experiments.marginal_lowdin_target import phase_alignment
from experiments.marginal_lowdin_target import build_target
import json, tempfile
from pathlib import Path

class TestLowdinTarget(unittest.TestCase):
    def test_phase_alignment(self):
        w=((1,0),(0,2)); g={w:F(2)}; s={w:F(-2)}
        self.assertEqual(phase_alignment(g,s,2)[2], (1,-1))
    def test_mismatch_is_visible(self):
        w=((1,0),(0,0)); self.assertGreater(phase_alignment({w:F(1)}, {w:F(2)}, 2)[0], 1e-7)
    def test_phase_helper_dimension_and_finite_inputs(self):
        with self.assertRaises(ValueError): phase_alignment({}, {}, 0)
        with self.assertRaises(ValueError): phase_alignment({}, {}, 9)

    def test_invalid_source_and_tolerance_refuse_before_pyscf(self):
        data={'modes':4,'particles':2,'geometry':[['H',[0,0,0]]],'basis':'sto-3g','unit':'Angstrom'}
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json'
            for changed in (dict(data,modes=True),dict(data,modes=4.0),dict(data,particles=True),dict(data,geometry=[])):
                source.write_text(json.dumps(changed))
                with self.assertRaises(ValueError): build_target(source,root/'out.json')
            source.write_text(json.dumps(data))
            for tol in (float('nan'),float('inf'),True,0,1e-3):
                with self.assertRaises(ValueError): build_target(source,root/'out.json',tolerance=tol)

if __name__ == '__main__': unittest.main()
