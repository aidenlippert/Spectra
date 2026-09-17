from fractions import Fraction as F
import copy
import math
import unittest

from research.intervention_reduction_20260916.exact import digest, sqrt_up
from research.intervention_reduction_20260916.test_exact import case
from research.intervention_reduction_20260916.trajectory import check_trajectory,polynomial_norm_range


class TrajectoryTest(unittest.TestCase):
    def euler(self):
        h=[[0,0,1,0],[0,0,0,0],[1,0,0,0],[0,0,0,0]]
        data,p=case(h,[[1,0],[0,0],[0,1]])
        t={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(p),
           'coefficient_denominator':10,'horizon_atomic_time':'1/10','phase_shift_Ha':'0',
           'observed_spatial_orbital':1,
           'segments':[{'duration_atomic_time':'1/10','controls_Ha':[],
                        'coefficients':[[[10,0],[0,0]],[[0,0],[0,-1]]]}]}
        return data,p,t

    def test_euler_against_independent_exact_two_level_evolution(self):
        data,p,t=self.euler();r=check_trajectory(data,p,t)
        # Full literal amplitudes are cos(t)|0> - i sin(t)|2>.
        lo,hi=map(float,map(F,r['population_interval']))
        self.assertLessEqual(lo,math.sin(.1)**2)
        self.assertGreaterEqual(hi,math.sin(.1)**2)
        self.assertEqual(F(r['segments'][0]['integrated_residual_squared']),F(1,30000))
        self.assertEqual(F(r['predicted_population']),F(1,101))
        self.assertEqual(F(r['jump_total']),0)

    def test_segment_jumps_are_charged(self):
        data,p,t=self.euler();t['segments']*=2;t['horizon_atomic_time']='1/5'
        r=check_trajectory(data,p,t)
        self.assertEqual(F(r['jump_total']),F(1,10))
        self.assertGreater(F(r['normalized_state_error_bound']),F(1,10))

    def test_phase_shift_and_joint_cancellation(self):
        data,p=case([[1000000,0],[0,1000000]],[[3]])
        t={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(p),
           'coefficient_denominator':1,'horizon_atomic_time':'100','phase_shift_Ha':'1000000',
           'observed_spatial_orbital':0,
           'segments':[{'duration_atomic_time':'100','controls_Ha':[],'coefficients':[[[1,0]]]}]}
        r=check_trajectory(data,p,t)
        self.assertEqual(F(r['normalized_state_error_bound']),0)
        self.assertEqual(r['population_interval'],['1','1'])
        self.assertEqual(F(r['uniform_in_time_state_error_bound']),0)

    def test_norm_enclosure_includes_interior_not_just_endpoints(self):
        # p(theta)=1-4 theta+4 theta^2 vanishes inside although endpoints are 1.
        lo,hi=polynomial_norm_range([[1]],[[[1,0]],[[-4,0]],[[4,0]]],1,1)
        self.assertLessEqual(lo,0)
        self.assertGreaterEqual(hi,1)
        data,p=case([[0,0],[0,0]],[[1]])
        t={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(p),
           'coefficient_denominator':1,'horizon_atomic_time':'1','phase_shift_Ha':'0',
           'observed_spatial_orbital':0,'segments':[{'duration_atomic_time':'1',
           'controls_Ha':[],'coefficients':[[[1,0]],[[-4,0]],[[4,0]]]}]}
        r=check_trajectory(data,p,t)
        self.assertIsNone(r['uniform_in_time_state_error_bound'])
        self.assertFalse(r['uniform_in_time_target_met'])

    def test_time_binding_and_malformed_coefficients_refuse(self):
        data,p,t=self.euler();bad=copy.deepcopy(t);bad['horizon_atomic_time']='1/5'
        with self.assertRaisesRegex(ValueError,'horizon'):
            check_trajectory(data,p,bad)
        bad=copy.deepcopy(t);bad['segments'][0]['coefficients'][0][0][0]=10.0
        with self.assertRaisesRegex(ValueError,'integer'):
            check_trajectory(data,p,bad)
        bad=copy.deepcopy(t);bad['proposal_sha256']='wrong'
        with self.assertRaisesRegex(ValueError,'binding'):
            check_trajectory(data,p,bad)

    def test_control_outside_box_refuses(self):
        h=[[0,1],[1,0]]
        data,p=case(h,[[1]],[(h,F(1,10))])
        t={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(p),
           'coefficient_denominator':1,'horizon_atomic_time':'1','phase_shift_Ha':'0',
           'observed_spatial_orbital':0,
           'segments':[{'duration_atomic_time':'1','controls_Ha':['1/5'],'coefficients':[[[1,0]]]}]}
        with self.assertRaisesRegex(ValueError,'amplitude box'):
            check_trajectory(data,p,t)


if __name__=='__main__':
    unittest.main()
