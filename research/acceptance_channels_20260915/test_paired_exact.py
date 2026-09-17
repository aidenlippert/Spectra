import copy
from fractions import Fraction
from itertools import combinations
import unittest
from experiments.marginal_symbolic import expand_squares as original
from research.acceptance_channels_20260915.paired_exact import expand_squares,CALLS
from research.acceptance_channels_20260915.test_channel import act,compose


def blocks():
    words=[((1,2*p),(0,2*q+1),(0,2*r+1)) for p in range(4) for q,r in combinations(range(4),2)]
    rows=[[((i*7+3)%11)-5 for i in range(len(words))],[(i%5)-2 for i in range(len(words))]]
    return [{'words':words,'factor':rows},
        {'words':[tuple((1-c,k) for c,k in reversed(w)) for w in words],'factor':copy.deepcopy(rows)}]


class ExactPairedTest(unittest.TestCase):
    def test_all_256_local_fock_columns_against_independent_ladder_action(self):
        sample=blocks();den=7;polynomial,_=expand_squares(sample,den,8)
        labels=list(range(8))
        operators=[]
        for row in sample[0]['factor']:
            B={w:Fraction(c,den) for w,c in zip(sample[0]['words'],row) if c}
            adjoint={tuple((1-c,k) for c,k in reversed(w)):v for w,v in B.items()}
            operators.append((B,adjoint))
        for state in range(256):
            expected={}
            for B,dagger in operators:
                for left,right in ((B,dagger),(dagger,B)):
                    for target,value in compose(left,right,state,labels).items():
                        expected[target]=expected.get(target,Fraction(0))+value
            self.assertEqual(act(polynomial,state,labels),{k:v for k,v in expected.items() if v})

    def test_dense_pair_matches_original_polynomial_and_statistics_exactly(self):
        sample=blocks();expected=original(sample,10**9,8)
        self.assertEqual(expand_squares(sample,10**9,8),expected)
        self.assertEqual(CALLS[-1]['exact_matched_block_pairs'],1)
        self.assertTrue(all(len(w)<=4 for w in expected[0]))

    def test_different_adjoint_coefficients_fall_back_without_assuming_pairing(self):
        sample=blocks();sample[1]['factor'][0][0]+=1
        self.assertEqual(expand_squares(sample,17,8),original(sample,17,8))
        self.assertEqual(CALLS[-1]['exact_matched_block_pairs'],0)

    def test_integer_products_have_no_machine_overflow(self):
        sample=blocks()
        for b in sample:b['factor']=[[v*10**40 for v in b['factor'][0]]]
        self.assertEqual(expand_squares(sample,10**41,8),original(sample,10**41,8))

    def test_duplicate_words_and_reversed_annihilators_reproduce_original(self):
        sample=blocks()
        w=sample[0]['words'][0];sample[0]['words'].append((w[0],w[2],w[1]))
        sample[1]['words'].append(tuple((1-c,k) for c,k in reversed(sample[0]['words'][-1])))
        for b in sample:
            for row in b['factor']:row.append(3)
        self.assertEqual(expand_squares(sample,19,8),original(sample,19,8))

    def test_original_refusal_paths_remain(self):
        for mutation in ('denominator','coefficient','word','length','charge'):
            sample=blocks();den=13
            if mutation=='denominator':den=True
            elif mutation=='coefficient':sample[1]['factor'][0][0]=.5
            elif mutation=='word':sample[1]['words'][0]=((1,100),(1,3),(0,0))
            elif mutation=='length':sample[1]['factor'][0].pop()
            else:sample[1]['words'][0]=((0,0),)
            with self.assertRaises(ValueError):expand_squares(sample,den,8)
            with self.assertRaises(ValueError):original(sample,den,8)


if __name__=='__main__':unittest.main()
