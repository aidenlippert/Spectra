from fractions import Fraction as F
import copy
import unittest
from research.intervention_reduction_20260916.kernel import CheckedKernel
from research.intervention_reduction_20260916.trajectory import check_trajectory
from research.intervention_reduction_20260916 import test_trajectory
from research.intervention_reduction_20260916.test_exact import case
from research.intervention_reduction_20260916.exact import digest


class KernelTest(unittest.TestCase):
    def test_compiled_residual_equals_literal_action_with_jumps_and_shift(self):
        data,p,t=test_trajectory.TrajectoryTest().euler()
        t['segments']*=2;t['horizon_atomic_time']='1/5';t['phase_shift_Ha']='7/3'
        direct=check_trajectory(data,p,t)
        kernel=CheckedKernel(data,p,1);compiled=check_trajectory(None,None,t,_kernel=kernel)
        ignored={'exact_replay_seconds','peak_RSS_bytes'}
        self.assertEqual({k:v for k,v in direct.items() if k not in ignored},
                         {k:v for k,v in compiled.items() if k not in ignored})
        self.assertEqual(kernel.receipt['configuration_records_retained_for_queries'],0)
        bad=copy.deepcopy(t);bad['observed_spatial_orbital']=0
        with self.assertRaisesRegex(ValueError,'kernel binding'):
            check_trajectory(data,p,bad,_kernel=kernel)

    def test_noncommuting_controls_nonorthogonal_basis_and_rational_weights(self):
        h=[[3,2,0,0],[2,1,4,0],[0,4,5,0],[0,0,0,2]]
        d=[[1,0,0,0],[0,-1,0,0],[0,0,2,0],[0,0,0,-2]]
        w=[[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]
        data,p=case(h,[[2,1],[1,-1],[0,1]],[(d,F(1)),(w,F(1))]);p['denominator']=7
        t={'kind':'rational_molecular_trajectory_v1','proposal_sha256':digest(p),
           'coefficient_denominator':100,'horizon_atomic_time':'1/10','phase_shift_Ha':'7/3',
           'observed_spatial_orbital':1,'segments':[{'duration_atomic_time':'1/10',
           'controls_Ha':['1/7','-2/3'],'coefficients':[[[100,0],[0,0]],[[1,-1],[2,1]],[[0,1],[1,-2]]]}]}
        direct=check_trajectory(data,p,t);kernel=CheckedKernel(data,p,1)
        compact=check_trajectory(None,None,t,_kernel=kernel)
        for key in ('segments','normalized_state_error_bound','uniform_in_time_state_error_bound','population_interval'):
            self.assertEqual(direct[key],compact[key])


if __name__=='__main__':
    unittest.main()
