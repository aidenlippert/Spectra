from fractions import Fraction
import unittest
from experiments.marginal_symbolic import add, adj, mono, product, scale
from research.acceptance_channels_20260915.augment import paired_entry
from research.acceptance_channels_20260915.channel import supplied_channel


class PairedBlockTest(unittest.TestCase):
    def test_all_cross_terms_match_a_joint_positive_square(self):
        B, _, _ = supplied_channel()
        words = list(B)
        for coefficients in ([2,-3,5,-7],[-1,2,0,4]):
            operator = add(*(mono(w,Fraction(v)) for w,v in zip(words,coefficients)))
            expected = add(product(adj(operator),operator),product(operator,adj(operator)))
            actual = add(*(scale(paired_entry(w,v),Fraction(a*b))
                for w,a in zip(words,coefficients) for v,b in zip(words,coefficients)))
            self.assertEqual(actual,expected)
            self.assertLessEqual(max(map(len,actual)),4)


if __name__ == '__main__':
    unittest.main()
