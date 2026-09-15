import unittest
from fractions import Fraction as F

from research.certificate_scaling.dressed_quadratic_recognition import dress_polynomial, recognize
from research.certificate_scaling.test_dressed_quadratic_certificate import occupation_matrix
from experiments.marginal_symbolic import canonical


def model(m, matching):
    q={():F(7)}
    for i in range(m): q[((1,i),(0,i))]=F(i+1)
    hubs=(m-2,m-1)
    edges=list(zip(range(0,m-2,2),range(1,m-2,2)))
    edges += [(h,v) for h in hubs for v in range(m-2)]
    for i,j in edges:
            q[((1,i),(0,j))]=q[((1,j),(0,i))]=F(1)
    return q, dress_polynomial(q, matching, m)


class RecognitionTests(unittest.TestCase):
    def test_m4_exact_and_pure_quadratic(self):
        q,p=model(4,((0,1),))
        r=recognize(p,4)
        self.assertTrue(r["accepted"]); self.assertTrue(r["exact_regeneration_match"])
        self.assertTrue(recognize(q,4)["accepted"])

    def test_permuted_two_hubs_sizes(self):
        for m in (4,6,8):
            q,p=model(m,tuple((i,i+1) for i in range(0,m-2,2)))
            perm=list(range(m)); perm[0],perm[2]=perm[2],perm[0]
            pp=dress_polynomial({tuple((a,perm[i]) for a,i in w):c for w,c in p.items()},(),m)
            r=recognize(pp,m)
            self.assertTrue(r["accepted"], (m,r))

    def test_refusals(self):
        _,p=model(4,((0,1),))
        bad=dict(p); w=next(w for w in bad if len(w)==4); bad[w]+=1
        self.assertFalse(recognize(bad,4)["accepted"])
        density=dict(p); density[((1,0),(1,1),(0,0),(0,1))]=F(1)
        self.assertFalse(recognize(density,4)["accepted"])
        self.assertFalse(recognize({((1,0),(0,1)):F(1)},4)["accepted"])
        self.assertFalse(recognize({((1,0),(0,1)):1.0},4)["accepted"])

    def test_independent_occupation_matrix_CZ(self):
        q,h=model(4,((0,1),))
        states,Q=occupation_matrix(q,4)
        _,H=occupation_matrix(h,4)
        signs=[-1 if s&1 and s&2 else 1 for s in states]
        self.assertEqual(H,[[signs[i]*Q[i][j]*signs[j] for j in range(16)] for i in range(16)])

    def test_disconnected_and_noncanonical_refusal(self):
        q={((1,0),(0,2)):F(1),((1,2),(0,0)):F(1)}
        h=dress_polynomial(q,((0,1),),4)
        self.assertEqual(recognize(h,4)['reason'],'vertex-deletion-disconnected')
        self.assertFalse(recognize({((0,0),(1,0)):F(1)},4)['accepted'])
        for coeff in (True,1.0,'1'):
            with self.assertRaises(ValueError): dress_polynomial({():coeff},(),4)
        with self.assertRaises(ValueError): dress_polynomial({((2,0),):F(1)},(),4)

    def test_invalid_matching(self):
        with self.assertRaises(ValueError): dress_polynomial({},((0,1),(1,2)),4)
        with self.assertRaises(ValueError): dress_polynomial({},((0,4),),4)


if __name__ == '__main__': unittest.main()
