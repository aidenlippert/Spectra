import unittest
from fractions import Fraction as F
from experiments.v7_certificate import Generator,check_certificate
from experiments.v7_adaptive_taylor import adaptive_taylor
class AdaptiveTaylorTests(unittest.TestCase):
 def test_direct_proof_checked_independently(self):
  for gamma in (F(0),F(1,5),F(2)):
   h={'ZI':F(1,2),'XX':F(1,3)};o={'YI':F(1)};T=F(1,5);tol=F(1,1000)
   p,w,c=adaptive_taylor(Generator(h,gamma,2),o,T,tol)
   r=check_certificate(Generator(h,gamma,2),o,[p],w,tol,expected_time=T)
   self.assertEqual(r['status'],'certified')
 def test_zero_generator_and_failure(self):
  p,w,c=adaptive_taylor(Generator({},F(0),1),{'X':F(1)},F(1),F(0))
  self.assertEqual(check_certificate(Generator({},F(0),1),{'X':F(1)},[p],w,F(0),expected_time=F(1))['status'],'certified')
  with self.assertRaises(ValueError):adaptive_taylor(Generator({'Z':F(1)},F(0),1),{'X':F(1)},F(1),F(0),max_order=2)
if __name__=='__main__':unittest.main()
