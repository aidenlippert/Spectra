import unittest
from fractions import Fraction as F
import numpy as np
from experiments.v6_headroom import Trace
from experiments.v6_audit import excitation,local_dependencies


class AuditTests(unittest.TestCase):
    def test_constant_input_obstruction_without_looking_at_outputs(self):
        train=[Trace(np.tile([50.,0.],(80,1)),np.zeros(81))]
        test=[Trace(np.tile([60.,0.],(80,1)),np.full(81,np.nan))]
        r=excitation(train,test)
        self.assertEqual(r['training_input_contrast_rank'],0)
        self.assertEqual(r['blocks_with_unidentified_input_contrast'],r['planned_blocks'])
        # Distinct input gains are exactly indistinguishable at the training input.
        u0=F(50);b=F(1,10);d=F(2);change=F(3,100)
        self.assertEqual(b*u0+d,(b+change)*u0+d-change*u0)
        self.assertNotEqual(b*60+d,(b+change)*60+d-change*u0)

    def test_forced_cayley_hamilton_coefficients_exactly(self):
        A=((F(1,2),F(1,5)),(F(1,5),F(1,3)))
        B=(F(1),F(2));x=(F(2),F(-1));u0=F(3);u1=F(-2)
        def step(state,u):return tuple(sum(A[i][j]*state[j] for j in range(2))+B[i]*u for i in range(2))
        x1=step(x,u0);x2=step(x1,u1)
        c1=F(-5,6);c0=F(19,150);D0=F(1,15);D1=F(1)
        self.assertEqual(x2[0]+c1*x1[0]+c0*x[0],D0*u0+D1*u1)

    def test_v4_calibration_is_not_a_v5_import_dependency(self):
        graph=local_dependencies(('experiments.v5_run','experiments.v5_account'))
        self.assertIn('experiments.v5_core',graph)
        self.assertTrue(all(not name.startswith('experiments.v4') for name in graph))


if __name__=='__main__':unittest.main()
