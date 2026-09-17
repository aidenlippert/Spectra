from fractions import Fraction
from itertools import combinations
import unittest
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from experiments.marginal_symbolic import add,mono,scale
from research.acceptance_channels_20260915.augment import paired_entry
from research.acceptance_channels_20260915.dense_t2 import tensors,contracted_entry,normal_word,embed_parent_gram
from research.acceptance_channels_20260915.rank_update import RankUpdatedNormal,factorize


class DenseT2Test(unittest.TestCase):
    def test_cross_contractions_against_full_literal_CAR(self):
        words=[((1,2*p),(0,2*q+1),(0,2*r+1)) for p in range(3) for q,r in combinations(range(3),2)]
        left=np.array([2,-1,3,0,1,-2,4,1,-3]);right=np.array([-1,2,1,3,-2,1,0,-1,2])
        expected=add(*(scale(paired_entry(w,v),Fraction(int(a*b))) for w,a in zip(words,left) for v,b in zip(words,right)))
        actual={}
        for word,value in contracted_entry(tensors(words,left,3),tensors(words,right,3),3):
            w,sign=normal_word(word)
            if w is not None:actual[w]=actual.get(w,Fraction(0))+Fraction(sign*value,2)
        actual={w:v for w,v in actual.items() if v}
        self.assertEqual(actual,expected)

    def test_normal_update_and_inverse_match_explicit_spd_matrix(self):
        rng=np.random.default_rng(23)
        base=sparse.diags(np.arange(1.,9.));C=sparse.csc_matrix(rng.normal(size=(8,3)))
        normal=RankUpdatedNormal(base,C)+sparse.eye(8)*.01
        dense=base.toarray()+C.toarray()@C.toarray().T+np.eye(8)*.01
        rhs=rng.normal(size=8)
        np.testing.assert_allclose(normal@rhs,dense@rhs,atol=1e-12)
        np.testing.assert_allclose(factorize(normal,splu).solve(rhs),np.linalg.solve(dense,rhs),atol=1e-12)

    def test_unsafe_integer_contraction_refused(self):
        with self.assertRaises(ValueError):tensors([((1,0),(0,1),(0,3))],[10**10],2)

    def test_repeated_words_checked_before_fixed_width_accumulation(self):
        word=((1,0),(0,1),(0,3))
        with self.assertRaises(ValueError):tensors([word,word],[2**62,2**62],2)
        C,A=tensors([word,word],[2**70,-2**70+3],2)
        self.assertEqual(C[0,0],3)

    def test_truncated_or_noninteger_coefficients_refused(self):
        word=((1,0),(0,1),(0,3))
        with self.assertRaises(ValueError):tensors([word],[],2)
        with self.assertRaises(ValueError):tensors([word],[.5],2)

    def test_nested_extension_preserves_the_entire_old_square(self):
        basis=np.array([[1,0,3],[0,2,1],[2,1,0]])
        parent={'basis_integers':basis[:,:2].tolist()}
        old=np.array([[2.,.3],[.3,1.]])
        new=embed_parent_gram(basis,parent,old,3)
        np.testing.assert_allclose(basis@new@basis.T,basis[:,:2]@old@basis[:,:2].T)
        self.assertGreaterEqual(np.linalg.eigvalsh(new).min(),0)
        with self.assertRaises(ValueError):embed_parent_gram(-basis,parent,old,3)


if __name__=='__main__':unittest.main()
