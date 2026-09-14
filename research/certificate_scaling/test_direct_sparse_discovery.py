"""Independent controls for the H-only sparse discovery experiment.

These tests check physical reference inequalities and rejection paths.  They
do not reimplement the LP or its CAR maps.
"""
from __future__ import annotations

import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np

from experiments.marginal_coefficient import hopping_model
from experiments.marginal_symbolic import mono
from research.certificate_scaling.direct_sparse_discovery import discover
from research.certificate_scaling.direct_pair_sdp import run as pair_run


class DirectSparseControls(unittest.TestCase):
    def test_constant_shift_and_hopping_known_ground(self):
        # Four modes, N=2, one matched hopping pair.  The one-particle bond
        # gives a known many-body ground energy -1; adding 3I shifts it to 2.
        hop = {
            ((1, 0), (0, 1)): F(-1),
            ((1, 1), (0, 0)): F(-1),
        }
        shifted = dict(hop)
        shifted[()] = F(3)
        cert, receipt = discover(shifted, 4, 2, budget=32, denominator=10**7)
        self.assertLessEqual(float(F(receipt["lower"])), 2.0 + 1e-12)
        self.assertEqual(receipt["source_factors_used"], False)
        self.assertEqual(receipt["source_upper_used"], False)

    def test_fixed_number_operator_has_known_floor(self):
        h = {((1, i), (0, i)): F(1) for i in range(4)}
        _, receipt = discover(h, 4, 2, budget=32, denominator=10**7)
        # Every N=2 determinant has energy exactly 2.
        self.assertLessEqual(float(F(receipt["lower"])), 2.0 + 1e-12)

    def test_quantized_sos_can_leave_nonzero_physical_residual(self):
        h = hopping_model(4, F(1, 5), False)
        _, receipt = discover(h, 4, 2, budget=32, denominator=10**6)
        self.assertGreater(F(receipt["residual_l1"]), 0)
        self.assertGreater(F(receipt["lower"]), F(receipt["numeric_proposed_b"]) - F(1, 1000))

    def test_ideal_body_two_replays_generated_degree_six(self):
        h = hopping_model(6, F(1, 5), False)
        _, receipt = discover(h, 6, 3, budget=64, ideal_body=2, denominator=10**6)
        self.assertGreaterEqual(receipt["maximum_coefficient_degree"], 6)
        self.assertEqual(receipt["source_factors_used"], False)
        self.assertEqual(receipt["source_upper_used"], False)
        self.assertLessEqual(float(F(receipt["lower"])), float(receipt["numeric_proposed_b"]) + 1e-8)

    def test_rejects_nonhermitian_hamiltonian(self):
        h = {((1, 0), (0, 1)): F(1)}
        with self.assertRaises(ValueError):
            discover(h, 4, 2, budget=16)

    def test_pair_sdp_matches_or_improves_same_budget_lower(self):
        # Both methods receive the same H-derived candidate budget.  The SDP
        # cone contains the rank-one pair choices, so its numerical proposal
        # should not be weaker; compare the exact replayed lower endpoints.
        h = hopping_model(4, F(1, 5), False)
        with tempfile.TemporaryDirectory() as td:
            _, lp = discover(h, 4, 2, budget=32, denominator=10**6)
            sdp = pair_run(h, 4, 2, 32, Path(td), denominator=10**6)
        # Rounded factors and separate numerical solvers can move the exact
        # replayed endpoints slightly; this run differs by <1e-6 Ha.
        self.assertGreaterEqual(
            float(F(sdp["lower"])) - float(F(lp["lower"])), -1e-5
        )


if __name__ == "__main__":
    unittest.main()
