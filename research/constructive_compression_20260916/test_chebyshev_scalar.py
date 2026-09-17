"""Exact diagonal examples exercise the gate and corrupted residual claims."""
from fractions import Fraction as F
import unittest
from research.constructive_compression_20260916.chebyshev_scalar import gate


def exact_vectors(energies,seed,b,ell,residuals):
    a=[(b+ell-2*e)/(b-ell) for e in energies]
    old=list(map(F,seed));previous=None
    for j,r in enumerate(residuals):
        new=[a[i]*old[i]+r[i] if j==0 else 2*a[i]*old[i]-previous[i]+r[i]
             for i in range(len(a))]
        previous,old=old,new
    return old


class TestGate(unittest.TestCase):
    def test_certifies_known_ground_with_an_upper_norm(self):
        v=exact_vectors([F(0),F(2)],[1,1],F(2),F(0),[[F(0),F(0)]]*8)
        out=gate(2,0,F(1,2),sum(map(abs,v)),[0]*8)
        self.assertTrue(out['accepted']);self.assertLessEqual(out['lower'],0)

    def test_projection_erasure_is_caught_by_residual_charge(self):
        # H=diag(-4,2), seed=(1,0), A=diag(5,-1); force every v_j=0.
        # Then r1=(-5,0), r2=(1,0), and all later residuals vanish.
        residuals=[[F(-5),F(0)],[F(1),F(0)]]+[[F(0),F(0)]]*6
        v=exact_vectors([F(-4),F(2)],[1,0],F(2),F(0),residuals)
        self.assertEqual(v,[0,0])
        self.assertFalse(gate(2,0,F(1,2),0,[5,1]+[0]*6)['accepted'])
        # Omitting errors really would give a false claim: data must be replayed.
        self.assertTrue(gate(2,0,F(1,2),0,[0]*8)['accepted'])

    def test_overlap_and_strict_boundary(self):
        self.assertFalse(gate(2,0,F(1,2),1,[0],overlap=1)['accepted'])
        self.assertTrue(gate(2,0,F(1,2),F(1,10),[0],overlap=F(1,2))['accepted'])

    def test_bad_premises_and_unchecked_floats(self):
        with self.assertRaises(TypeError):gate(2,0,.5,0,[0])
        with self.assertRaises(ValueError):gate(1,2,F(1,2),0,[0])
        with self.assertRaises(ValueError):gate(2,0,F(1,2),-1,[0])
        with self.assertRaises(ValueError):gate(2,0,F(1,2),0,[-1])
        with self.assertRaises(ValueError):gate(2,0,F(1,2),0,[])
        with self.assertRaises(ValueError):gate(2,0,F(1,2),0,[0],overlap=0)
