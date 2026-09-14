"""Small exact controls for the rational independence check."""
from fractions import Fraction as F
import unittest
from research.mechanism_transfer_20260913.ledger import modular_rank


class RationalRank(unittest.TestCase):
    def test_row_swap_fractions_and_known_linear_dependence(self):
        a=[F(0),F(1,3),F(2,3)];b=[F(1,2),F(0),F(5,2)]
        # The first two columns have determinant -1/6.
        self.assertEqual(modular_rank([a,b]),2)
        self.assertEqual(modular_rank([a,b,[2*x-3*y for x,y in zip(a,b)]]),2)

    def test_noninvertible_denominator_is_refused(self):
        with self.assertRaises(ValueError):modular_rank([[F(1,2**61-1)]])


if __name__=='__main__':unittest.main()
