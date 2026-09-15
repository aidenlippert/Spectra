"""Exact replay checks; no numerical optimizer is invoked by these tests."""
from fractions import Fraction as F
import unittest

import numpy as np
from experiments.marginal_hopping import (
    SECTOR, hamiltonian, seed_blocks, verify_certificate,
    baseline_blocks, solve, export_certificate,
)


def exact_seed_certificate():
    trial = [int(s == 11) for s in SECTOR]  # sites 0,1,3 occupied
    return {"mode_count": 6, "particle_number": 3, "denominator": 1,
            "b_numerator": 1, "trial": trial,
            "blocks": [{"name": b["name"], "words": b["words"], "factor": [[1]]}
                       for b in seed_blocks()]}


class HoppingTests(unittest.TestCase):
    def test_exact_zero_hopping_certificate(self):
        result = verify_certificate(hamiltonian(0), exact_seed_certificate())
        self.assertEqual(result["lower"], "1")
        self.assertEqual(result["upper"], "1")
        self.assertEqual(result["residual_row_norm"], "0")

    def test_wrong_scalar_is_charged_to_residual(self):
        c = exact_seed_certificate()
        c["b_numerator"] = 1001
        result = verify_certificate(hamiltonian(0), c)
        self.assertEqual(result["lower"], "1")
        self.assertEqual(result["residual_row_norm"], "1000")

    def test_hopping_residual_is_not_discarded(self):
        result = verify_certificate(hamiltonian(F(1, 5)), exact_seed_certificate())
        self.assertEqual(result["lower"], "2/5")
        self.assertEqual(result["upper"], "1")
        self.assertEqual(result["width"], "3/5")

    def test_corrupted_square_changes_safe_bound(self):
        c = exact_seed_certificate()
        c["blocks"][0]["factor"] = [[2]]
        result = verify_certificate(hamiltonian(0), c)
        self.assertGreater(F(result["residual_row_norm"]), 0)
        self.assertLessEqual(F(result["lower"]), 1)

    def test_malformed_certificate_rejected(self):
        for key, value in (("denominator", 0), ("particle_number", 2),
                           ("trial", [0]*20), ("b_numerator", 1.5)):
            c = exact_seed_certificate()
            c[key] = value
            with self.assertRaises(ValueError):
                verify_certificate(hamiltonian(0), c)
        c = exact_seed_certificate()
        c["blocks"][0]["factor"] = [[0.5]]
        with self.assertRaises(ValueError):
            verify_certificate(hamiltonian(0), c)
        with self.assertRaises(ValueError):
            verify_certificate(np.asarray(hamiltonian(0), dtype=float), exact_seed_certificate())

    def test_matched_spectrum_formula(self):
        for t in (0, .05, .2, 1, 2):
            spectrum = np.linalg.eigvalsh(np.asarray(hamiltonian(F(str(t))), dtype=float))
            actual = spectrum[0]
            expected = 2-t-np.sqrt(1+2*t+4*t*t)
            self.assertAlmostEqual(actual, expected, places=12)
            predicted = [1-t]*8 + [1+t]*8
            predicted += [2-t+s*np.sqrt(1+2*t+4*t*t) for s in (-1, 1)]
            predicted += [2+t+s*np.sqrt(1-2*t+4*t*t) for s in (-1, 1)]
            np.testing.assert_allclose(spectrum, sorted(predicted), atol=1e-12, rtol=1e-12)

    def test_solver_export_matches_known_quadratic_bound(self):
        h = hamiltonian(0)
        blocks = baseline_blocks()
        proposal = solve(h, blocks)
        certificate, _ = export_certificate(h, blocks, proposal)
        result = verify_certificate(h, certificate)
        self.assertGreater(F(result["lower"]), F(749, 1000))
        self.assertLessEqual(F(result["lower"]), F(3, 4))
        self.assertEqual(F(result["upper"]), 1)


if __name__ == "__main__":
    unittest.main()
