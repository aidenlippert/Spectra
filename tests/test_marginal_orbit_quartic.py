import unittest
import numpy as np
from fractions import Fraction as F

from experiments.marginal_orbit_quartic import (symmetry_group, word_action,
    run, average_solution)
from experiments.marginal_coefficient import dictionaries, hopping_model, solve_coefficients


class OrbitQuarticTests(unittest.TestCase):
    def test_group_and_action(self):
        g = symmetry_group(6)
        self.assertEqual(len(g), 12)
        self.assertTrue(all(word_action(((0, 0),), p)[0] in
                            [((0, i),) for i in range(6)] for p in g))

    def test_projection_is_psd(self):
        b = dictionaries(4, "mixed")
        s = solve_coefficients(hopping_model(4, F(1, 5)), 4, 2, b)
        a = average_solution(s, b, 4)
        self.assertTrue(all(np.linalg.eigvalsh(q).min() > -1e-7 for q in a["grams"]))

    def test_exact_replay(self):
        c, r = run(4)
        self.assertLess(float(r["width_float"]), 1.0)
        self.assertIn("residual_l1", r)


if __name__ == "__main__":
    unittest.main()
