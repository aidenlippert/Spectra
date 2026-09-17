from fractions import Fraction as F
import unittest
from research.constructive_response_20260916.bond_coercivity import (
    bond_b_upper,graph_bound,verify_bond_all_sectors,bond_matrix,
)
from research.constructive_response_20260916.composition_exact import is_psd


class BondTests(unittest.TestCase):
    def test_all_occupation_landmarks(self):
        self.assertEqual(bond_b_upper(F(0)),2)
        self.assertEqual(bond_b_upper(F(3)),1)
        self.assertEqual(bond_b_upper(F(8)),1)

    def test_exact_PSD_in_all_16_states(self):
        for t in (F(0),F(1,2),F(1),F(2)):
            for a in (0,1,2,3,4,8):
                ok,b=verify_bond_all_sectors(F(a),t)
                self.assertTrue(ok,(a,t,b))

    def test_understated_allowance_rejected(self):
        self.assertFalse(is_psd(bond_matrix(F(0),F(1),F(199,100))))
        self.assertFalse(is_psd(bond_matrix(F(3),F(1),F(99,100))))
        # D is a count, not the bit values 2 and 8.
        M=bond_matrix(F(1),F(0),F(0))
        self.assertEqual(M[3][3],1)
        self.assertEqual(M[12][12],1)
        self.assertEqual(M[15][15],2)

    def test_overlap_and_refusal(self):
        c,b=graph_bound(F(8),F(1),4,F(1))
        self.assertEqual(c,4);self.assertGreater(b,1)
        with self.assertRaises(ValueError):
            verify_bond_all_sectors(F(-1))


if __name__=='__main__':
    unittest.main()
