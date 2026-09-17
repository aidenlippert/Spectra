from fractions import Fraction as F
import unittest
from research.constructive_compression_20260916.tensor_exact import inner,apply_mpo,sqrt_upper,recurrence_residual


def v(a,b):return {'widths':[1,1],'denominators':[1],'layers':[[[0,0,0,a],[0,1,0,b]]]}
H={'widths':[1,1],'layers':[[[0,0,[-4,0,0,2],'1']]]}

class ExactTensorTests(unittest.TestCase):
    def test_exact_noninteger_scaling_and_mpo(self):
        a=v(3,4);a['denominators']=[5]
        self.assertEqual(inner(a,a),1)
        self.assertEqual(inner(a,apply_mpo(a,H)),F(-4,25))
        self.assertEqual(inner(apply_mpo(a,H),apply_mpo(a,H)),F(208,25))

    def test_residual_matches_known_vector_and_detects_erasure(self):
        n,e=recurrence_residual(v(0,0),v(1,0),None,H,2,0)
        self.assertEqual(n,25);self.assertEqual(e,5)
        n,e=recurrence_residual(v(0,0),v(0,0),v(1,0),H,2,0)
        self.assertEqual(n,1)
        n,e=recurrence_residual(v(5,0),v(1,0),None,H,2,0)
        self.assertEqual(n,0)

    def test_sqrt_is_an_upper_bound_exactly(self):
        for x in (F(0),F(2),F(1,7),F(10**40,3),F(1,10**50)):
            u=sqrt_upper(x,40);self.assertGreaterEqual(u*u,x)
            if u:self.assertLess((u-F(1,1<<40))**2,x)
