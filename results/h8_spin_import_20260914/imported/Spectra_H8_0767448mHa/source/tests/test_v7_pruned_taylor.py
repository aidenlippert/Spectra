import unittest
from fractions import Fraction as F
from experiments.v7_certificate import Generator,check_certificate
from experiments.v7_pruned_taylor import pruned_taylor
class PrunedTaylorTests(unittest.TestCase):
 def test_pruning_proof_includes_omissions(self):
  h={'ZI':F(1,2),'XX':F(1,10000),'IZ':F(1,3)};o={'YI':F(1)};T=F(1,5);tol=F(1,1000)
  for mode in ('magnitude','work'):
   p,w,c=pruned_taylor(Generator(h,F(1,5),2),o,T,tol,mode=mode)
   self.assertGreater(F(c['pruning_spent']),0)
   r=check_certificate(Generator(h,F(1,5),2),o,[p],w,tol,expected_time=T)
   self.assertEqual(r['status'],'certified')
 def test_omitted_bound_cannot_be_removed(self):
  h={'ZI':F(1,2),'XX':F(1,10000)};o={'YI':F(1)};T=F(1,5);tol=F(1,1000)
  p,w,c=pruned_taylor(Generator(h,F(0),2),o,T,tol)
  for key,val in w['witnesses'].items():
   if key!='jump:0' and val:w['witnesses'][key]=[];break
  self.assertEqual(check_certificate(Generator(h,F(0),2),o,[p],w,tol,expected_time=T)['status'],'rejected')
if __name__=='__main__':unittest.main()
