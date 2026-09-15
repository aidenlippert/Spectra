import unittest

from experiments.stabilizer import diagnose_css


class StabilizerDiagnosticTests(unittest.TestCase):
    def test_bridge_outside_buffer_is_witness(self):
        # X1X2 and X2X3 generate X1X3, supported in A={1,3},
        # but neither generator is confined to B=A.
        x, z = diagnose_css([(1, 1, 0), (0, 1, 1)], [], 3, {0, 2}, {0, 2})
        self.assertFalse(x.passed)
        self.assertEqual(x.witness, (1, 0, 1))
        self.assertIsNone(z.witness)

    def test_buffer_containing_bridge_passes(self):
        x, _ = diagnose_css([(1, 1, 0), (0, 1, 1)], [], 3, {0, 2}, {0, 1, 2})
        self.assertTrue(x.passed)

    def test_small_repetition_like_patch_passes(self):
        x, z = diagnose_css([(1, 1, 0), (0, 1, 1)], [], 3, {1}, {0, 1, 2})
        self.assertTrue(x.passed)
        self.assertTrue(z.passed)

    def test_rejects_noncommuting_css_generators(self):
        with self.assertRaises(ValueError):
            diagnose_css([(1, 0)], [(1, 0)], 2, {0}, {0})

    def test_rejects_bad_region_and_dimensions(self):
        with self.assertRaises(ValueError):
            diagnose_css([(1, 0)], [], 2, {1}, {0})
        with self.assertRaises(ValueError):
            diagnose_css([(1, 2)], [], 2, {0}, {0})
        with self.assertRaises(ValueError):
            diagnose_css([(1, 0)], [], 2, {0.5}, {0.5})
        with self.assertRaises(ValueError):
            diagnose_css([(1, 0)], [], 2.0, {0}, {0})


if __name__ == "__main__":
    unittest.main()
