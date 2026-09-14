import unittest
from experiments.scaling import bounded_closure, quadratic_generators, xx_label, run

class ScalingTests(unittest.TestCase):
    def test_complete_result_has_verified_zero_residual(self):
        rows = run(max_n=2)
        self.assertTrue(all(row['verified_zero_residual'] for row in rows))

    def test_xx_labels(self):
        self.assertEqual(xx_label(4, 1), "IXXI")

    def test_small_complete_and_cap(self):
        controls = quadratic_generators(2)
        seeds = [next(iter(h)) for h in controls]
        basis, complete, _, _ = bounded_closure(controls, seeds, cap=64)
        self.assertTrue(complete)
        self.assertGreaterEqual(len(basis), len(seeds))
        _, complete, _, omissions = bounded_closure(controls, seeds, cap=1)
        self.assertFalse(complete)
        self.assertTrue(omissions)
        with self.assertRaises(ValueError): bounded_closure(controls, seeds, cap=0)
        with self.assertRaises(ValueError): bounded_closure(controls, seeds, work_cap=-1)

if __name__ == "__main__": unittest.main()
