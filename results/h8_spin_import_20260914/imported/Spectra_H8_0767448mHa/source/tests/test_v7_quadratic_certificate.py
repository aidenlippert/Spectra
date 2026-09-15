import unittest
from copy import deepcopy
from fractions import Fraction as F
from experiments.v7_certificate import Generator, Piece
from experiments.v7_quadratic_certificate import check,adaptive_with_square
from experiments.v7_headroom import model,TOL
class QuadraticTests(unittest.TestCase):
 def test_nontrivial_tighter_norm(self):
  # R=XI+IX+ZZ, R²=3I+2XX, norm=sqrt(5), better than sqrt(2)+1.
  o={'XI':F(1),'IX':F(1),'ZZ':F(1)}
  pieces=[Piece(F(1),(o,))]
  from experiments.v7_quadratic_probe import square
  from experiments.v7_certificate import norm_witness
  from experiments.certificates import _sqrt_interval
  s,_=square(o);self.assertEqual(s,{'II':F(3),'XX':F(2)})
  gs,_=norm_witness(s);hi=_sqrt_interval(F(5),16)[1]
  w=dict(schema='v7-quadratic-1',integration_basis='power',claimed_bound=str(hi),witnesses={'jump:0':dict(kind='square',groups=gs,upper=str(hi)),'residual:0:0':dict(kind='groups',groups=[])})
  self.assertEqual(check(Generator({},F(0),2),{},pieces,w,F(23,10),expected_time=F(1))['status'],'certified')
  bad=deepcopy(w);bad['witnesses']['jump:0']['upper']='2'
  self.assertEqual(check(Generator({},F(0),2),{},pieces,bad,F(3),expected_time=F(1))['status'],'rejected')
 def test_development_degree_saving_is_checked(self):
  h,o=model(3,'xxz');T=F(1,2)
  p,w,c=adaptive_with_square(Generator(h,F(2),3,max_terms=512),o,T,TOL)
  self.assertEqual(len(p.coefficients)-1,11)
  self.assertTrue(any(v['kind']=='square' for v in w['witnesses'].values()))
  self.assertEqual(check(Generator(h,F(2),3,max_terms=512),o,[p],w,TOL,expected_time=T)['status'],'certified')
 def test_partial_square_partition_and_corruption(self):
  h,o=model(3,'xxz');T=F(1,2)
  p,w,c=adaptive_with_square(Generator(h,F(2),3,max_terms=512),o,T,TOL,trigger_ratio=F(11,10),block_groups=4)
  self.assertEqual(len(p.coefficients)-1,11)
  self.assertEqual(check(Generator(h,F(2),3,max_terms=512),o,[p],w,TOL,expected_time=T)['status'],'certified')
  part=next(v for v in w['witnesses'].values() if v['kind']=='sum')
  part['parts'][0]['labels'].append(part['parts'][0]['labels'][0])
  self.assertEqual(check(Generator(h,F(2),3,max_terms=512),o,[p],w,TOL,expected_time=T)['status'],'rejected')
if __name__=='__main__':unittest.main()
