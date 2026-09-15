from fractions import Fraction as F
import unittest
from experiments.v2_chain import predicted_risk, enumerate_risk, source_eigenvalues


class ChainTests(unittest.TestCase):
    def test_full_menu_risks_and_increasing_marginals(self):
        for r in range(1, 5):
            risks = [enumerate_risk(r, k) for k in range(r + 1)]
            self.assertEqual(risks, [predicted_risk(r, k) for k in range(r + 1)])
            gains = [risks[k] - risks[k + 1] for k in range(r)]
            self.assertTrue(all(b == 2 * a for a, b in zip(gains, gains[1:])))

    def test_density_spectrum_and_zero_visibility_boundary(self):
        for r in range(1, 7):
            for v in (F(0), F(1, 2), F(1)):
                eig = source_eigenvalues(r, v)
                self.assertEqual(sum(value * mult for value, mult in eig), 1)
                self.assertTrue(all(value >= 0 for value, _ in eig))
            self.assertEqual(enumerate_risk(r, r, F(0)), F(1, 2))

    def test_resource_and_type_bounds(self):
        for args in ((7, 0), (0, 0), (3, 4), (True, 0), (3, 1.5)):
            with self.assertRaises(ValueError): enumerate_risk(*args)


if __name__ == '__main__': unittest.main()
