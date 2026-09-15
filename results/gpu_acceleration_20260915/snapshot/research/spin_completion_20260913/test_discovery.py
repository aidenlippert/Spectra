"""Cross-check numerical maps against exact CAR, and frame inclusion."""
from fractions import Fraction as F
from pathlib import Path
import json
import unittest
import numpy as np

from research.certificate_scaling.commutator_dual_witness import moment_decode
from research.joint_patterns_20260913.dual import integer_grams
from research.spin_subspace_20260913.full_dual import frames as old_frames
from research.spin_completion_20260913.discovery import Model


class CompletionMaps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        prior = root/'results/molecular_collective_20260913/campaign/h6'
        cls.model = Model(json.loads((prior/'fixture.json').read_text()), json.loads((prior/'rank_10/tail.json').read_text()))
        cls.dual = json.loads((root/'results/spin_subspace_20260913/full_dual/witness.json').read_text())

    def test_full_frame_contains_every_previous_operator(self):
        m = self.model
        _, old = old_frames(m.p, m.tail, m.symmetry['parity_masks'])
        all_polys = [p for g in m.frames for p in g['polynomials']]
        self.assertEqual(len(all_polys), 492)
        self.assertEqual([len(g['polynomials']) for g in m.frames], [30, 30, 93, 93, 93, 93, 30, 30])
        for group in old:
            for p in group['polynomials']:
                self.assertIn(p, all_polys)
        for C in m.frame_coefficients:
            self.assertEqual(np.linalg.matrix_rank(C, tol=1e-10), C.shape[1])

    def test_new_charge_and_mixed_charge_maps_match_exact_moments(self):
        m = self.model; y = moment_decode(self.dual['moments'], 12)
        values = np.array([float(y.get(w, F(0))) for w in m.rows])
        for gid in (0, 2):
            gram, den = integer_grams(m.frames[gid]['polynomials'], [y], 'anticommutator')
            exact = np.array([[float(F(c, den)) for c in row] for row in gram[0]])
            C = m.frame_coefficients[gid]
            actual = C.T@(m.price_maps[gid]@values).reshape(C.shape[0], C.shape[0])@C
            self.assertTrue(np.allclose(actual, exact, atol=2e-14, rtol=2e-14))


if __name__ == '__main__':
    unittest.main()
