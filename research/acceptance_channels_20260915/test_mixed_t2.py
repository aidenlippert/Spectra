import unittest
from itertools import combinations
from fractions import Fraction
import numpy as np
from experiments.marginal_symbolic import add,scale
from research.acceptance_channels_20260915.augment import paired_entry
from research.acceptance_channels_20260915.dense_t2 import normal_word
from research.acceptance_channels_20260915.mixed_t2 import general_tensor,general_entry,operator_moment_matrix
from research.interacting_scaling_20260915.dictionary import representative
from research.interacting_scaling_20260915.singlet_trace import singlet_trace


class MixedT2Test(unittest.TestCase):
    def test_general_contractions_including_creator_annihilator_overlap(self):
        words=[((1,p),(0,q),(0,r)) for p in range(3) for q,r in combinations(range(3),2)]
        left=np.array([2,-1,3,0,1,-2,4,1,-3]);right=np.array([-1,2,1,3,-2,1,0,-1,2])
        expected=add(*(scale(paired_entry(w,v),Fraction(int(a*b))) for w,a in zip(words,left) for v,b in zip(words,right)))
        actual={}
        for word,value in general_entry(general_tensor(words,left,3),general_tensor(words,right,3),3):
            w,sign=normal_word(word)
            if w is not None:actual[w]=actual.get(w,Fraction(0))+Fraction(sign*value,2)
        self.assertEqual({w:v for w,v in actual.items() if v},expected)

    def test_vectorized_moments_match_exact_singlet_trace_of_each_pair(self):
        spin=lambda k:1 if k%2==0 else -1
        words=[((1,p),(0,q),(0,r)) for p in range(4) for q,r in combinations(range(4),2)
            if spin(p)-spin(q)-spin(r)==1]
        polynomials=[[paired_entry(a,b) for b in words] for a in words]
        keys={representative(((1,a),(0,b))) for a in range(4) for b in range(4) if spin(a)==spin(b)}
        keys.update(representative(((1,a),(1,b),(0,c),(0,d)))
            for a,b in combinations(range(4),2) for c,d in combinations(range(4),2)
            if spin(a)+spin(b)==spin(c)+spin(d))
        # Include zero-moment canonical words as supplied entries, rather than
        # silently filling a required unavailable value.
        moments={w:float(singlet_trace({w:Fraction(1)},4,2)) for w in keys}
        M=operator_moment_matrix(words,moments,4)
        expected=np.array([[float(singlet_trace(p,4,2)) for p in row] for row in polynomials])
        np.testing.assert_allclose(M,expected,atol=1e-14)

    def test_required_unsupplied_moment_is_refused(self):
        with self.assertRaises(ValueError):
            operator_moment_matrix([((1,0),(0,0),(0,1))],{},4)


if __name__=='__main__':unittest.main()
