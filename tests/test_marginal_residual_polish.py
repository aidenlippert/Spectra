import unittest
from fractions import Fraction as F
from experiments.marginal_residual_polish import polish_block
from experiments.marginal_signed_atoms import verify_block, residual


class ResidualPolishTests(unittest.TestCase):
    def test_fixed_nonconstant_metric_pair_rescaling_recovers_exact_boundary(self):
        # A=[1,4][1,4]^T-2I; the new residual scale differs from W.
        a=[[F(-1),F(4)],[F(4),F(14)]]
        initial={'metric_weights':[1,2],'atom_scale':1,'atoms':[]}
        result,bound,d=polish_block(a,initial)
        self.assertEqual(bound,-2);self.assertTrue(d['accepted'])
        self.assertEqual(d['pair_edges_processed'],1)
        verify_block(a,result,-2,True)
        with self.assertRaises(ValueError): verify_block(a,result,F(-199,100),True)

    def test_disconnected_residual_components_get_independent_positive_scales(self):
        a=[[F(0) for _ in range(4)] for _ in range(4)]
        a[0][0]=-1;a[0][1]=a[1][0]=4;a[1][1]=14
        a[2][2]=4;a[3][3]=9;a[2][3]=a[3][2]=-6
        initial={'metric_weights':[1]*4,'atom_scale':1,'atoms':[]}
        result,bound,d=polish_block(a,initial)
        self.assertAlmostEqual(float(bound),-2,places=7)
        self.assertTrue(all(w>0 for w in d['scaling_weights']))
        self.assertEqual(d['pair_edges_processed'],2);verify_block(a,result,bound,True)

    def test_existing_pair_direction_is_merged_without_duplicate_atoms(self):
        a=[[F(1),F(4)],[F(4),F(16)]]
        initial={'metric_weights':[1,1],'atom_scale':2,
                 'atoms':[{'indices':[0,1],'amplitudes':[1,4],'amplitude_scale':4,'weight':16}]}
        result,bound,d=polish_block(a,initial)
        self.assertEqual(bound,0);self.assertEqual(len(result['atoms']),1)
        self.assertEqual(residual(a,result['metric_weights'],result['atoms'],result['atom_scale'],True)[0],[[0,0],[0,0]])


if __name__=='__main__': unittest.main()
