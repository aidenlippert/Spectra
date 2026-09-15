import copy
from fractions import Fraction as F
from itertools import product as tuples
import unittest
import numpy as np
from experiments.marginal_symbolic import word_product,expand_squares
from experiments.marginal_coefficient import dagger,gram_map,coefficient_rows
from research.correlated_pair_20260913.test_mps import correlated_fixture,explicit_moment
from research.reconstruction_compression_20260914.moments import Proposal
from research.reconstruction_compression_20260914.reduced import project,cut_maps
from research.reconstruction_compression_20260914.ideal_rank import rank_mod,certify
from experiments.marginal_symbolic import multiplier_basis

class ReducedTests(unittest.TestCase):
    def test_proposal_moments_with_fermionic_signs(self):
        data,state=correlated_fixture();s=Proposal(data,state)
        letters=list(tuples((0,1),range(4)))
        words=[w for d in (0,1,2,3,4) for w in tuples(letters,repeat=d)]
        s.fill(words)
        for w in words:self.assertAlmostEqual(s.cache[w],float(explicit_moment(w)),places=13)

    def test_congruence_all_degrees_against_exact_factors(self):
        words=[((1,0),(0,2),(0,1)),((1,3),(0,1),(0,0)),((0,0),)]
        rows=coefficient_rows(4);lookup={w:i for i,w in enumerate(rows)}
        V=np.array([[1,2],[3,-1],[2,0]],dtype=float);R=np.array([[2,1],[1,-1]])
        Q=R.T@R;P,ii,jj=project(gram_map(words,lookup),V)
        blocks=[{'words':words,'factor':(R@V.T).astype(int).tolist()}]
        exact,_=expand_squares(blocks,1,4,max_degree=3)
        self.assertTrue(any(len(w)==6 and c for w,c in exact.items()))
        np.testing.assert_array_equal(P@Q[ii,jj],np.array([float(exact.get(w,0)) for w in rows]))

    def test_odd_adjoint_pair_cancels_sextic_terms(self):
        words=[((1,0),(0,2),(0,1)),((1,3),(0,1),(0,0))]
        blocks=[{'words':ws,'factor':[[2,3]]} for ws in (words,list(map(dagger,words)))]
        poly,_=expand_squares(blocks,1,4,max_degree=3)
        self.assertFalse(any(len(w)==6 for w in poly))

    def test_cut_congruence_keeps_physical_PSD(self):
        words=[((1,k),(0,j),(0,i)) for i in range(4) for j in range(i+1,4) for k in range(4)]
        rng=np.random.default_rng(43);A=rng.normal(size=(24,len(words)));C=A.T@A
        V=cut_maps(words,C,4)
        np.testing.assert_allclose(V.T@V,np.eye(4),atol=1e-12)
        B=rng.normal(size=(5,len(words)));M=B.T@B
        self.assertGreaterEqual(np.linalg.eigvalsh(V.T@M@V)[0],-1e-10)

    def test_modular_injectivity_and_refusal(self):
        self.assertEqual(rank_mod([{0:2,1:3},{0:4,1:6}],2)[0],1)
        self.assertEqual(rank_mod([{0:2,1:3},{0:4,1:7}],2)[0],2)
        with self.assertRaises(ValueError):rank_mod([{0:1}],1,9)
        proof=certify(6,3,multiplier_basis(6,max_body=2))
        self.assertEqual(proof['rank'],proof['columns'])
        # At four modes the exterior multiplication can have a kernel; this
        # fixture ensures the reducer refuses rather than assuming a theorem.
        with self.assertRaises(ValueError):certify(4,2,multiplier_basis(4,max_body=2))

    def test_proposal_does_not_bypass_structure_or_zero_norm(self):
        data,base=correlated_fixture()
        for what in ('binding','charge','zero','duplicate'):
            state=copy.deepcopy(base)
            if what=='binding':state['fixture_sha256']='bad'
            if what=='charge':state['tensors'][0][0][1]=1
            if what=='zero':state['tensors'][1]=[]
            if what=='duplicate':state['tensors'][0].append(state['tensors'][0][0])
            with self.subTest(what=what),self.assertRaises(ValueError):Proposal(data,state)

if __name__=='__main__':unittest.main()
