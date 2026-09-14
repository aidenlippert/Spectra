from fractions import Fraction as F
import unittest
import numpy as np
from experiments.marginal_symbolic import mono,add,scale,canonical,adj,product
from research.certificate_scaling.spin_twirl import twirl
from research.collective_completion_20260914.spin_share import share


def flip(p):
    return canonical({tuple((c,i^1) for c,i in w):v for w,v in p.items()})


def square(p):
    return product(adj(p),p)


class SpinShareTests(unittest.TestCase):
    def test_spin_flip_invariance_with_fermion_signs(self):
        p=add(mono(((1,0),(0,2),(0,3)),F(2,3)),
              mono(((1,2),(0,0),(0,1)),F(-7,5)))
        self.assertEqual(twirl(square(p)),twirl(square(flip(p))))

    def test_self_orbit_cross_parity_vanishes(self):
        words=[((1,i),(0,j)) for i in range(4) for j in range(4) if i%2==j%2]
        groups=[{'name':'spin conserving quadratic','words':words}]
        maps,ids,receipt=share(groups,[([(0,np.eye(len(words)))],'quadratic')])
        self.assertEqual(ids,[0,0])
        self.assertEqual(sum(m[0][0][1].shape[1] for m in maps),len(words))
        factors=[]
        for members,name in maps:
            B=members[0][1]
            factors.append([add(*(scale(mono(w),F(str(float(v)))) for w,v in zip(words,col))) for col in B.T])
        for a in factors[0]:
            for b in factors[1]:
                self.assertEqual(twirl(add(product(adj(a),b),product(adj(b),a))),{})

    def test_distinct_orbit_preserves_averaged_positive_terms(self):
        words=[((0,0),),((0,2),)];other=[((0,1),),((0,3),)]
        groups=[{'name':'alpha','words':words},{'name':'beta','words':other}]
        maps,ids,receipt=share(groups,[([(0,np.eye(2))],'alpha'),([(1,np.eye(2))],'beta')])
        self.assertEqual(ids,[0]);np.testing.assert_array_equal(maps[0][0][0][1],np.eye(2))
        p=add(mono(words[0],F(2)),mono(words[1],F(-3)))
        q=add(mono(other[0],F(5)),mono(other[1],F(7)))
        self.assertEqual(twirl(add(square(p),square(q))),twirl(add(square(p),square(flip(q)))))


if __name__=='__main__':unittest.main()
