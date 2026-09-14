import unittest
from fractions import Fraction as F
from experiments.v11_guarded_taylor import guarded_taylor,norm_rejection
from experiments.v8_integer_taylor import fraction_free_taylor
from experiments.v7_certificate import Generator,check_certificate
from experiments.v7_headroom import model,TOL


def ledger():return dict(rejected_orders=0,inspected_max_entries=0,squared_entries=0,max_comparison_bits=0)


class GuardTests(unittest.TestCase):
    def test_strict_boundary_zero_and_distinct_guards(self):
        # Weight1, derivative coefficient1, tolerance1: equality is not a refusal.
        for mode in ('max','frobenius','cascade'):
            self.assertIsNone(norm_rejection({'X':1},1,F(1),F(1),0,mode,ledger()))
            self.assertIsNone(norm_rejection({},1,F(1),F(0),0,mode,ledger()))
            self.assertTrue(norm_rejection({'X':1},1,F(1),F(0),0,mode,ledger()))
        self.assertIsNone(norm_rejection({'X':1,'Y':1},1,F(1),F(1),0,'max',ledger()))
        self.assertEqual(norm_rejection({'X':1,'Y':1},1,F(1),F(1),0,'frobenius',ledger()),'normalized_Hilbert_Schmidt')

    def test_integer_comparison_equals_rational_inequality(self):
        for k in range(5):
            for tol in (F(0),F(1,1000),F(1,7),F(1),F(100)):
                for time in (F(1,5),F(3,7),F(2)):
                    nums={'XI':-7,'YY':3,'ZZ':11};denom=17
                    bound_square=(time**(k+1)/F(k+1))**2*sum(F(v,denom)**2 for v in nums.values())
                    gate=norm_rejection(nums,denom,time,tol,k,'frobenius',ledger())
                    self.assertEqual(bool(gate),bound_square>tol**2)

    def test_identical_outputs_and_independent_checks(self):
        for n,fam,gamma,time in [(3,'xxz',F(2),F(1,2)),(3,'mixed',F(1,5),F(1,5)),(4,'xxz',F(0),F(1,5))]:
            h,o=model(n,fam)
            p0,w0,_=fraction_free_taylor(Generator(h,gamma,n,512),o,time,TOL)
            for mode in ('max','frobenius','cascade'):
                p,w,c=guarded_taylor(Generator(h,gamma,n,512),o,time,TOL,guard=mode)
                self.assertEqual(p,p0);self.assertEqual(w['witnesses'],w0['witnesses'])
                self.assertEqual(w['claimed_bound'],w0['claimed_bound'])
                self.assertGreater(c['guard_work']['rejected_orders'],0)
                self.assertEqual(check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=time)['status'],'certified')

    def test_stationary_and_invalid_mode(self):
        p,w,_=guarded_taylor(Generator({},0,1),{'Z':F(1,3)},F(1),F(0))
        self.assertEqual(p.coefficients,({'Z':F(1,3)},));self.assertEqual(w['claimed_bound'],'0')
        with self.assertRaises(ValueError):guarded_taylor(Generator({},0,1),{'Z':F(1)},F(1),F(0),guard='unknown')


if __name__=='__main__':unittest.main()
