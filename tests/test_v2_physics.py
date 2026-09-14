import unittest
from fractions import Fraction
import numpy as np
from experiments.v2_physics import World, born_probability, density

class PhysicsTests(unittest.TestCase):
    def test_density_and_born(self):
        for a in (0,1):
            for b in (0,1):
                for s in (-1,1):
                    rho=density(a,b,s); self.assertAlmostEqual(np.trace(rho).real,1)
                    self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(),-1e-12)
                    for action in ((0,0),(0,1),(1,0),(1,1)):
                        i,j=action; pa=np.array([[0,1],[1,0]]) if i==0 else np.diag([1,-1]); pb=np.array([[0,1],[1,0]]) if j==0 else np.diag([1,-1])
                        self.assertAlmostEqual(np.trace((np.eye(4)+np.kron(pa,pb))/2@rho).real,float(born_probability(a,b,s,action)))
    def test_parity_shared_by_calibration_and_target(self):
        w=World((1,0,1),(0,1,1),Fraction(1),Fraction(0)); rng=np.random.default_rng(1); ca=(1,0,1); cb=(0,1,0)
        self.assertEqual(set(w.calibration("a",ca,8,rng)),{0}); self.assertEqual(set(w.calibration("b",cb,8,rng)),{1})
        self.assertTrue(all(x==1 for x in w.target(ca,cb,1,(0,1),8,rng)))
    def test_validation(self):
        with self.assertRaises(ValueError): World((1,),(0,1))
        with self.assertRaises(ValueError): density(0,0,1,Fraction(3,2))
        with self.assertRaises(ValueError): World((1,),(0,),Fraction(1),Fraction(0)).target((1,0),(0,),1,(0,0),1,np.random.default_rng())
if __name__=="__main__": unittest.main()
