import unittest
from fractions import Fraction as F
from itertools import combinations
from experiments.marginal_symbolic import canonical,mono,add
from experiments.marginal_spin_reduction import spin_operators,commutator
from research.certificate_scaling.spin_basis import decompose_words,spin_action,spin_weight2

class SpinBasisTests(unittest.TestCase):
    def test_local_spin_action_matches_independent_full_commutator(self):
        p=canonical(add(mono(((1,1),(0,0),(0,3)),F(3,7)),mono(((1,2),(0,4),(0,5)),F(-2,3))))
        sp,sm,_=spin_operators(6)
        self.assertEqual(spin_action(p),commutator(sp,p))
        self.assertEqual(spin_action(p,False),commutator(sm,p))

    def test_complete_mixed_cubic_decomposition_and_character_counts(self):
        for m in (4,6):
            words=[((0,i),) for i in range(m)]+[((1,k),(0,i),(0,j)) for i,j in combinations(range(m),2) for k in range(m)]
            groups,r=decompose_words(words,m)
            counts={j:sum(len(g['copies']) for g in groups if g['two_spin']==j) for j in (1,3)}
            d1=sum(spin_weight2(w)==1 for w in words);d3=sum(spin_weight2(w)==3 for w in words)
            self.assertEqual(counts,{1:d1-d3,3:d3});self.assertEqual(r['descendant_dimension'],len(words))
            self.assertLessEqual(r['max_local_dimension'],8)
            for g in groups:
                for ladder in g['copies']:
                    self.assertFalse(spin_action(ladder[0]))

    def test_nonclosed_span_and_noncommuting_parity_refuse(self):
        with self.assertRaises(ValueError):decompose_words([((1,0),)],4)
        with self.assertRaises(ValueError):decompose_words([((1,0),),((1,1),)],4,[1])
        groups,r=decompose_words([((1,i),) for i in range(4)],4,[3])
        self.assertEqual(len(groups),2);self.assertEqual(r['invariant_Gram_entries'],2)

if __name__=='__main__':unittest.main()
