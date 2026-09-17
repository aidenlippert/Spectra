from fractions import Fraction as F
import math
import unittest
from experiments.marginal_symbolic import encode,mono,add
from research.intervention_reduction_20260916.exact import digest
from research.intervention_reduction_20260916.test_exact import case
from research.intervention_reduction_20260916.robust import control_norm_bound,check_robust


class ControlNeighborhoodTest(unittest.TestCase):
    def test_nonadjacent_local_control_norm(self):
        hopping=add(mono(((1,0),(0,5)),1),mono(((1,5),(0,0)),1))
        bound,meta=control_norm_bound(encode(hopping),6)
        self.assertEqual(bound,1)
        self.assertEqual(meta['support_modes'],2)
        self.assertEqual(meta['local_dimension'],4)

    def test_constant_perturbation_against_analytic_evolution(self):
        h=[[0]*4 for _ in range(4)];w=[[0]*4 for _ in range(4)];w[0][2]=w[2][0]=1
        data,p=case(h,[[1]],[(w,F(1,10))])
        t={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(p),
           'coefficient_denominator':1,'horizon_atomic_time':'1','phase_shift_Ha':'0',
           'observed_spatial_orbital':1,'segments':[{'duration_atomic_time':'1',
           'controls_Ha':['0'],'coefficients':[[[1,0]]]}]}
        r=check_robust(data,p,t,[F(1,100)])
        self.assertEqual(F(r['nominal_state_error_bound']),0)
        self.assertEqual(F(r['normalized_state_error_bound']),F(1,100))
        lo,hi=r['population_interval_float']
        self.assertLessEqual(lo,math.sin(.01)**2)
        self.assertGreaterEqual(hi,math.sin(.01)**2)
        with self.assertRaisesRegex(ValueError,'nonnegative'):
            check_robust(data,p,t,[F(-1,100)])
        with self.assertRaisesRegex(ValueError,'per control'):
            check_robust(data,p,t,[])


if __name__=='__main__':
    unittest.main()
