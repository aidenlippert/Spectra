import unittest
from fractions import Fraction as F
from experiments.v7_certificate import Generator, derive_certificate, check_certificate
from experiments.v7_reducers import full_taylor, projected_taylor, bfs_basis, residual_basis, arnoldi

class ReducerTests(unittest.TestCase):
 def test_rotation_polynomials(self):
  make=lambda:Generator({'Z':F(1,2)},F(0),1)
  constructors=[lambda g:full_taylor(g,{'X':F(1)},F(1,10),6),lambda g:projected_taylor(g,{'X':F(1)},F(1,10),6,{'X','Y'}),lambda g:residual_basis(g,{'X':F(1)},F(1,10),6,cap=2),lambda g:arnoldi(g,{'X':F(1)},F(1,10),6,2)]
  for constructor in constructors:
   piece,cost=constructor(make())
   self.assertAlmostEqual(float(piece.coefficients[1]['Y']),-1)
   self.assertAlmostEqual(float(piece.coefficients[2]['X']),-.5)
   w=derive_certificate(make(),{'X':F(1)},[piece])
   r=check_certificate(make(),{'X':F(1)},[piece],w,F(1,10**8),expected_time=F(1,10))
   self.assertEqual(r['status'],'certified')
 def test_basis_refusals_and_projection_error(self):
  g=lambda:Generator({'Z':F(1,2)},F(0),1)
  with self.assertRaises(RuntimeError): bfs_basis(g(),{'X':F(1)},2,cap=1)
  b,c=bfs_basis(g(),{'X':F(1)},2,cap=2)
  self.assertEqual(b,{'X','Y'})
  p,c=projected_taylor(g(),{'X':F(1)},F(1),4,{'X'})
  w=derive_certificate(g(),{'X':F(1)},[p])
  self.assertEqual(check_certificate(g(),{'X':F(1)},[p],w,F(1,10),expected_time=F(1))['status'],'over_tolerance')
if __name__=='__main__': unittest.main()
