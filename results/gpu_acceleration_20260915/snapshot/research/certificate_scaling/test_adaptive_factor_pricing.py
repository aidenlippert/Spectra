"""Independent controls for adaptive factor pricing."""
from __future__ import annotations

import json
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))

import numpy as np

from experiments.marginal_coefficient import hopping_model
from experiments.marginal_symbolic import verify
from research.certificate_scaling.adaptive_factor_pricing import pricing_vectors, run


class AdaptiveFactorPricingControls(unittest.TestCase):
    def test_pricing_finds_negative_sparse_direction(self):
        matrix = np.array([[0.0, 1.0 / 8.0], [1.0 / 8.0, 1.0 / 4.0]])
        proposals = pricing_vectors(matrix, width=2, count=4)
        self.assertTrue(proposals)
        self.assertLess(proposals[0][0], 0.0)
        self.assertEqual(sorted(proposals[0][1].tolist()), [0, 1])

    def test_small_run_replays_and_tracks_complete_map(self):
        h = hopping_model(4, F(1, 5), False)
        with tempfile.TemporaryDirectory() as td:
            cert, receipt = run(h, 4, 2, Path(td), width=2, iterations=2,
                                batch=8, atom_cap=512, ideal_body=1,
                                seconds=30, denominator=10**6)
            checked = verify(cert)
            history = json.loads((Path(td) / "history.json").read_text())
        self.assertEqual(checked["lower"], receipt["lower"])
        self.assertGreater(receipt["complete_pricing_map_nonzeros"], 0)
        self.assertEqual(receipt["source_factors_used"], False)
        self.assertEqual(receipt["source_upper_used"], False)
        self.assertFalse(receipt["omitted_family_optimality_proved"])
        self.assertEqual(receipt["rounds_solved"], len(history))
        self.assertGreaterEqual(receipt["coefficient_rows"], 1)

    def test_body_two_ideal_replays_higher_degree_residual(self):
        h = hopping_model(4, F(1, 5), False)
        with tempfile.TemporaryDirectory() as td:
            cert, receipt = run(h, 4, 2, Path(td), width=2, iterations=1,
                                batch=4, atom_cap=512, ideal_body=2,
                                seconds=30, denominator=10**6)
        self.assertEqual(verify(cert)["lower"], receipt["lower"])
        # Body-2 ideals require the complete degree-6 coefficient basis even
        # when the accepted atom support happens not to generate degree six.
        self.assertGreater(receipt["coefficient_rows"], 47)
        self.assertLessEqual(receipt["residual_max_degree"], 6)


if __name__ == "__main__":
    unittest.main()
