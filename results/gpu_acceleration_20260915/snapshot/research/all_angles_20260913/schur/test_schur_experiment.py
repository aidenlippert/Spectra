import unittest
from fractions import Fraction as F

try:
    from .schur_experiment import schur_lower, run
except ImportError:
    from schur_experiment import schur_lower, run


class SchurTest(unittest.TestCase):
    def test_endpoint_is_below_oracle(self):
        H = [[F(0), F(1, 10), F(0)],
             [F(1, 10), F(2), F(1, 5)],
             [F(0), F(1, 5), F(3)]]
        r = schur_lower(H, 0)
        self.assertEqual(r["endpoint"], F(-1, 180))
        self.assertLessEqual(float(r["endpoint"]), -0.0050208)

    def test_excited_eigenvector_refuses(self):
        r = schur_lower([[F(0), F(0)], [F(0), F(1)]], 1)
        self.assertEqual(r["status"], "refuse_mu_le_theta")
        self.assertIsNone(r["endpoint"])

    def test_h4_car_reconstruction_and_refusal(self):
        result = run()["molecular_H4_RHF_determinant"]
        self.assertEqual(result["dimension"], 70)
        self.assertEqual(result["status"], "refuse_mu_le_theta")
        self.assertLessEqual(float(F(result["mu"])), float(F(result["theta"])))
        self.assertEqual(result["full_complement_entries"], 69 * 69)


if __name__ == "__main__":
    unittest.main()
