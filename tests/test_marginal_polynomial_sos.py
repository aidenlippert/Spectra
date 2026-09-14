import copy
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from experiments.marginal_polynomial_sos import integer_psd, replay_separator
from experiments.marginal_symbolic import encode


class PolynomialSquareTests(unittest.TestCase):
    def test_rank_one_denominator_scaling_preserves_psd_and_refusals(self):
        for scale in (2,17,10**20+3):
            for base,weights in [([[2,1,0],[1,2,0],[0,0,0]],[2,-3,0]),
                                 ([[0,0,0],[0,1,0],[0,0,2]],[0,1,2]),
                                 ([[0,0,0],[0,0,0],[0,0,0]],[1,2,3])]:
                matrix=[[scale*x+weights[i]*weights[j] for j,x in enumerate(row)] for i,row in enumerate(base)]
                self.assertEqual(integer_psd(matrix,initial_divisor=scale),integer_psd(matrix))
            bad=[[scale,2*scale],[2*scale,scale]]
            with self.assertRaises(ValueError):integer_psd(bad,initial_divisor=scale)
        with self.assertRaises(ValueError):integer_psd([[1,0],[0,1]],initial_divisor=2)
        for divisor in (0,-1,True,F(1)):
            with self.assertRaises(ValueError):integer_psd([[1]],initial_divisor=divisor)

    def test_fraction_free_positive_singular_and_indefinite_matrices(self):
        cases=[([[0,0,0],[0,2,1],[0,1,2]],2),
               ([[2,2,2],[2,2,2],[2,2,3]],2),
               ([[1,2,3],[2,4,6],[3,6,9]],1),([[0]],0)]
        for matrix,rank in cases:
            result=integer_psd(matrix);self.assertEqual(result['rank'],rank)
            self.assertEqual(result['nullity'],len(matrix)-rank)
        for matrix in [[[-1]],[[0,1],[1,0]],[[1,2],[2,1]],
                       [[1,0],[1,1]],[[F(1)]],[[True]],[],[[1,2]]]:
            with self.assertRaises(ValueError):integer_psd(matrix)

    def test_complete_square_basis_and_actual_polynomial_binding(self):
        c={'modes':8,'particles':4,'hamiltonian':encode({}),
           'polynomial_metric':{'denominator':1,'terms':[{'powers':[0]*4,'coefficient':1}]},
           'kind':'joint_polynomial_sos_separator_v1','target_lower':'1','square_degree':2,
           'functional':{'states':[15],'weights':['1']}}
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No determinant actions')):
            result=replay_separator(c)
        self.assertEqual(result['functional_value'],'-1')
        self.assertEqual(result['moment_psd'],{'dimension':37,'rank':1,'nullity':36,'positive_semidefinite':True})
        c['charge_square_degree']=2
        localized=replay_separator(c)
        self.assertEqual(localized['charge_moment_psd']['rank'],1)
        self.assertEqual(localized['charge_moment_psd']['dimension'],37)
        for key,value in [('target_lower','-1'),('square_degree',True),('square_degree',4)]:
            bad=copy.deepcopy(c);bad[key]=value
            with self.assertRaises(ValueError):replay_separator(bad)
        for degree in [True,-1,3]:
            bad=copy.deepcopy(c);bad['charge_square_degree']=degree
            with self.assertRaises(ValueError):replay_separator(bad)


if __name__=='__main__':unittest.main()
