from fractions import Fraction as F
import unittest
from research.transfer_followup_20260915.family_membership import in_span


class Span(unittest.TestCase):
    def test_coupled_coefficients_are_not_independent_directions(self):
        a, b, c = ('a',), ('b',), ('c',)
        basis = [{a: F(2), b: F(2)}]
        self.assertTrue(in_span({a: F(1), b: F(1)}, basis))
        self.assertFalse(in_span({a: F(1)}, basis))
        self.assertFalse(in_span({c: F(1)}, basis))
        self.assertTrue(in_span({a: F(1)}, basis+[{b: F(3)}]))


if __name__ == '__main__':
    unittest.main()
