"""Independent controls for optional adaptive-pricing modes."""
from __future__ import annotations

import json
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_coefficient import hopping_model
from experiments.marginal_symbolic import hermitian, multiplier_basis
from research.certificate_scaling.adaptive_factor_pricing import run


class AdaptivePricingOptionsControls(unittest.TestCase):
    # For M4,t=1/5,N2 the symmetric trial amplitudes (1,3,1) have
    # norm20 and energy numerator -14/5: exact variational upper=-7/50.
    def test_full_initial_checkpoint_and_pruning_replay(self):
        h = hopping_model(4, F(1, 5), False)
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            _, receipt = run(h, 4, 2, out, iterations=1, batch=4,
                             atom_cap=512, ideal_body=1, seconds=30,
                             prune=True, half_rows=False, denominator=10**6)
            history = json.loads((out / "history.json").read_text())
        self.assertEqual(receipt["rounds_solved"], 2)
        self.assertEqual(receipt["prune_inactive_atoms"], True)
        self.assertGreater(receipt["cumulative_atoms_pruned"], 0)
        self.assertGreaterEqual(receipt["cumulative_atoms_added"],
                                receipt["atoms_retained"])
        self.assertLessEqual(float(receipt["lower_float"]), -0.14 + 1e-8)
        self.assertEqual(len(history), 2)

    def test_full_and_half_initial_objectives_and_replay(self):
        h = hopping_model(4, F(1, 5), False)
        values = []
        for half in (False, True):
            with tempfile.TemporaryDirectory() as td:
                out = Path(td)
                _, receipt = run(h, 4, 2, out, iterations=1, batch=4,
                                 atom_cap=512, ideal_body=1, seconds=30,
                                 half_rows=half, prune=False, denominator=10**6)
                history = json.loads((out / "history.json").read_text())
                values.append((history[0]["numeric_lower"], receipt))
        self.assertAlmostEqual(values[0][0], values[1][0], places=7)
        self.assertLessEqual(float(F(values[1][1]["lower"])), -0.14 + 1e-8)

    def test_multiplier_basis_is_hermitian_for_control_sector(self):
        basis = multiplier_basis(4, max_body=2)
        self.assertTrue(all(hermitian(q) for q in basis))


if __name__ == "__main__":
    unittest.main()
