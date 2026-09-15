import unittest
from fractions import Fraction as F
from selected_schur import inertia_psd, run

class ExactControls(unittest.TestCase):
    def test_indefinite_ldl_rejects(self):
        self.assertFalse(inertia_psd([[F(1),F(2)],[F(2),F(1)]], F(0)))
    def test_zero_coupling(self):
        self.assertTrue(inertia_psd([[F(2),F(0)],[F(0),F(3)]], F(2)))
    def test_original_fixture_and_support(self):
        r=run('results/certificate_scaling/active_space_ladder/h4/fixture.json',[15,30,45,51,60,75,90,102,105,120,135,150,153,165,180,195,204,210,225,240])
        self.assertEqual(r['coupling_frobenius_sq'], '0')

if __name__=='__main__': unittest.main()
