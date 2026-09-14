import unittest
from fractions import Fraction as F
from experiments.v7_certificate import Generator,check_certificate
from experiments.v7_headroom import model,TOL
from experiments.v7_reducers import projected_taylor as old_projected
from experiments.v9_galerkin import adaptive_galerkin,projected_taylor

class GalerkinTests(unittest.TestCase):
 def test_last_omission_and_exact_restart(self):
  g=Generator({'Z':F(1)},F(0),1,512)
  c,om=projected_taylor(g,{'X':F(1)},{'X'},F(1,2),0)
  self.assertEqual(om,({'Y':F(-2)},))
  cs,om=projected_taylor(g,{'X':F(1)},{'X','Y'},F(1,2),4)
  old,_=old_projected(Generator({'Z':F(1)},F(0),1,512),{'X':F(1)},F(1,2),4,{'X','Y'})
  self.assertEqual(cs,old.coefficients)
  self.assertEqual(om,({},)*5)

 def test_nontrivial_physical_certificates(self):
  cases=[(2,{'ZI':F(1,3),'XX':F(2,5),'IZ':F(1,7)},{'ZI':F(1)},F(1,5),F(1,5))]
  h,o=model(3,'xxz');cases.append((3,h,o,F(2),F(1,5)))
  for n,h,o,gamma,T in cases:
   result=adaptive_galerkin(Generator(h,gamma,n,512),o,T,TOL)
   self.assertEqual(result['status'],'certified',result)
   ans=check_certificate(Generator(h,gamma,n,512),o,[result['piece']],result['certificate'],TOL,expected_time=T)
   self.assertEqual(ans['status'],'certified')
   self.assertGreater(result['expansions'],0)
   ans=check_certificate(Generator(h,gamma,n,512),o,[result['piece']],result['certificate'],TOL,expected_time=T+1)
   self.assertEqual(ans['status'],'rejected')

 def test_refusal_and_request_validation(self):
  h,o=model(3,'mixed')
  ans=adaptive_galerkin(Generator(h,F(0),3,512),o,F(1,2),TOL,max_expansions=0)
  self.assertEqual(ans['status'],'refused');self.assertEqual(ans['reason'],'expansion budget')
  for kwargs in ({'max_expansions':17},{'max_order':25},{'batch':0}):
   with self.assertRaises(ValueError):adaptive_galerkin(Generator(h,F(0),3,512),o,F(1,2),TOL,**kwargs)
  with self.assertRaises(ValueError):projected_taylor(Generator(h,F(0),3,512),o,{'ABX'},F(1,2))

 def test_bfs_shared_engine(self):
  h,o=model(3,'xxz');T=F(1,5)
  r=adaptive_galerkin(Generator(h,F(2),3,512),o,T,TOL,growth='bfs')
  self.assertEqual(r['status'],'certified',r)
  replay=check_certificate(Generator(h,F(2),3,512),o,[r['piece']],r['certificate'],TOL,expected_time=T)
  self.assertEqual(replay['status'],'certified')
  with self.assertRaises(ValueError):adaptive_galerkin(Generator(h,F(2),3,512),o,T,TOL,growth='unknown')

if __name__=='__main__':unittest.main()
