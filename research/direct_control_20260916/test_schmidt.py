import unittest
from fractions import Fraction as F
import numpy as np
from .schmidt_obstruction import lower_sqrt,run
from .validated_tensor import right_canonicalize

class TestSchmidt(unittest.TestCase):
 def test_square_root_direction(self):
  for x in [F(0),F(1,3),F(1,100000),F(17,19)]:
   lo=lower_sqrt(x)
   self.assertLessEqual(lo*lo,x)
   self.assertGreaterEqual((lo+F(1,1<<79))**2,x)
 def test_known_schmidt_tail_with_complex_gauge(self):
  # |psi> = .8|00> + .6|11>, exact normalized rational reference.
  ts=[np.array([[[.8,0],[0,.6]]],complex),np.array([[[1],[0]],[[0],[1]]],complex)]
  canonical,e=right_canonicalize(ts)
  rec=run(canonical,e+1e-14,F(1),[1,2])
  low=F(rec['lower_bounds']['1']['state_norm_error_at_least'])
  self.assertLessEqual(low,F(3,5));self.assertGreater(float(low),.599999)
  self.assertEqual(F(rec['lower_bounds']['2']['state_norm_error_at_least']),0)

if __name__=='__main__':unittest.main()
