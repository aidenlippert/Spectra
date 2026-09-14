import unittest
from fractions import Fraction as F
from .lift import apply, lift, check

class LiftTests(unittest.TestCase):
    def test_tiny_independent_action(self):
        # One-body hop plus diagonal term; verify CAR signs and Q filtering.
        H={((1,2),(0,0)):F(2),((1,0),(0,0)):F(3),((1,0),(0,0)):F(3)}
        # Direct application from |0> gives 2|2> and 3|0> (duplicate key
        # intentionally tests the same sparse action convention).
        out=apply(H,{1:F(1)})
        self.assertEqual(out[1],F(3)); self.assertEqual(out[4],F(2))
    def test_frozen_certificates_have_sparse_support(self):
        import json
        from research.compact_response_20260913 import program
        for name,limit in (('h6',30),('h8',30)):
            data,_,_=program.load_case(name); c=lift(data, F(1,100)); self.assertEqual(check(data,c)['upper_Ha'],c['upper_Ha'])
            self.assertLessEqual(c['chi_terms'],limit)
            self.assertEqual(len(c['chi_support']),len(c['chi_coefficients']))
            self.assertGreaterEqual(F(c['denominator']),F(1))

    def test_check_rejects_fixture_or_moment_mutation(self):
        from copy import deepcopy
        from research.compact_response_20260913 import program
        data,_,_=program.load_case('h6'); c=lift(data,F(1,100)); bad=deepcopy(c);bad['alpha']='1/3'
        with self.assertRaises(ValueError): check(data,bad)

if __name__=='__main__': unittest.main()
