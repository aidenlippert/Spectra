from fractions import Fraction as F
import unittest
from experiments.marginal_asymmetric_reference import model
from experiments.marginal_symbolic import transform


class AsymmetricReferenceTests(unittest.TestCase):
    def test_only_actual_symmetries_are_used(self):
        h=model(F(1,1000));mapping=[(i+5)%10 for i in range(10)]
        self.assertEqual(h,transform(h,mapping))
        for word in h:
            charge=[0]*5
            for creation,i in word:charge[i%5]+=2*creation-1
            self.assertEqual(charge,[0]*5)
        swap=list(range(10));swap[0],swap[1]=swap[1],swap[0];swap[5],swap[6]=swap[6],swap[5]
        self.assertNotEqual(h,transform(h,swap))


if __name__=='__main__':unittest.main()
