"""Exact controls for tuple indexing, component minimum, and Gram replay."""
from fractions import Fraction as F
from itertools import combinations
import unittest
from research.certificate_scaling.wedge_spectral_bound import matrix_lower, body_components, replay_residual

class WedgeSpectralControls(unittest.TestCase):
    def test_gram_with_residual_and_indefinite_control(self):
        lower,_=matrix_lower([[2,1],[1,2]], {'ell':'1','factor':[[1],[1]],'denominator':1})
        self.assertEqual(lower,1)
        lower,_=matrix_lower([[0,1],[1,0]], {'ell':'-1','factor':[[1],[1]],'denominator':1})
        self.assertEqual(lower,-1)
        # Arbitrary factors cannot create a false positive lower endpoint.
        lower,_=matrix_lower([[0,1],[1,0]], {'ell':'10','factor':[[10],[10]],'denominator':1})
        self.assertLessEqual(lower,-1)

    def test_wedge_tuples_sign_and_dimension(self):
        residual={((1,0),(1,1),(0,0),(0,1)): F(-2),
                  ((1,0),(1,2),(0,0),(0,2)): F(-3),
                  ((1,0),(1,1),(0,0),(0,2)): F(-1),
                  ((1,0),(1,2),(0,0),(0,1)): F(-1)}
        parts=body_components(residual,4,2)
        self.assertEqual(sum(len(indices) for indices,_ in parts),6)
        joined=next((indices,matrix) for indices,matrix in parts if len(indices)==2)
        self.assertEqual(joined, ([(0,1),(0,2)], [[F(2),F(1)],[F(1),F(3)]]))

    def test_disconnected_blocks_use_minimum_not_sum(self):
        residual={((1,0),(0,0)):F(-1),((1,1),(0,1)):F(-2)}
        proofs=[]
        for indices,matrix in body_components(residual,2,1):
            proofs.append({'indices':[list(i) for i in indices],'ell':str(matrix[0][0]),'factor':[[]],'denominator':1})
        result=replay_residual({'modes':2,'particles':1,'b':'0'},residual,
            {'method':'fixed_N_wedge_spectral_gram_v1','bodies':[{'body':1,'components':proofs}]})
        self.assertEqual(F(result['lower']),-2)
        proofs[0]['indices']=[[1]]
        with self.assertRaises(ValueError):replay_residual({'modes':2,'particles':1,'b':'0'},residual,
            {'method':'fixed_N_wedge_spectral_gram_v1','bodies':[{'body':1,'components':proofs}]})

    def test_malformed_or_nonhermitian_rejected(self):
        for witness in ({'ell':'0','factor':[[1],[1]],'denominator':0},
                        {'ell':'0','factor':[[1],[1,2]],'denominator':1},
                        {'ell':'0','factor':[[True],[1]],'denominator':1}):
            with self.assertRaises(ValueError):matrix_lower([[1,0],[0,1]],witness)
        with self.assertRaises(ValueError):matrix_lower([[0,1],[0,0]],{'ell':'0','factor':[[],[]],'denominator':1})

    def test_cubic_sign_and_full_particle_lift(self):
        residual={tuple((1,i) for i in indices)+tuple((0,i) for i in indices):F(2)
                  for indices in combinations(range(4),3)}
        bodies=[]
        for k in (1,2,3):
            proofs=[]
            for indices,matrix in body_components(residual,4,k):
                self.assertEqual(matrix[0][0], -2 if k==3 else 0)
                proofs.append({'indices':[list(i) for i in indices],'ell':str(matrix[0][0]),'factor':[[]],'denominator':1})
            bodies.append({'body':k,'components':proofs})
        result=replay_residual({'modes':4,'particles':4,'b':'0'},residual,
            {'method':'fixed_N_wedge_spectral_gram_v1','bodies':bodies})
        self.assertEqual(F(result['lower']),-8)

if __name__=='__main__':unittest.main()
