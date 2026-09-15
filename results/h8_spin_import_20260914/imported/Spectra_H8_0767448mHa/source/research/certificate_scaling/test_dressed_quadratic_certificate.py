import copy
import unittest
from fractions import Fraction as F
from research.certificate_scaling.dressed_quadratic_recognition import dress_polynomial
from research.certificate_scaling.dressed_quadratic_certificate import propose, replay


def quadratic_fixture():
    q = {(): F(3)}
    for i, d in enumerate((-2,-1,2,3)):
        q[((1,i),(0,i))] = F(d)
    for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3)):
        q[((1,i),(0,j))] = q[((1,j),(0,i))] = F(1,4)
    return q


def occupation_matrix(h, modes, particles=None):
    """Independent right-to-left occupation-mask action, with fermion signs."""
    states = [s for s in range(1 << modes) if particles is None or s.bit_count() == particles]
    pos = {s:i for i,s in enumerate(states)}
    out = [[F(0) for _ in states] for _ in states]
    for col,s in enumerate(states):
        for word,c in h.items():
            t,sign = s,1
            for create,i in reversed(word):
                if bool(t & (1<<i)) == bool(create):
                    sign=0
                    break
                if (t & ((1<<i)-1)).bit_count()%2: sign=-sign
                t ^= 1<<i
            if sign and t in pos: out[pos[t]][col] += c*sign
    return states,out


class DressedQuadraticCertificateTests(unittest.TestCase):
    def test_round_trip_and_independent_sector(self):
        import numpy as np
        q=quadratic_fixture(); h=dress_polynomial(q,((0,1),),4)
        for n in (0,1,2,3,4):
            with self.subTest(particles=n):
                payload=propose(h,4,n,denominator=10**8)
                result=replay(payload)
                _,matrix=occupation_matrix(h,4,n)
                exact_numeric=np.linalg.eigvalsh(np.array(matrix,dtype=float))[0]
                self.assertLessEqual(result['lower_float'],exact_numeric+1e-12)
                self.assertGreaterEqual(result['upper_float'],exact_numeric-1e-12)
                self.assertLess(result['width_float'],1e-4)
                self.assertEqual(payload['matching'],[[0,1]])

    def test_refuses_perturbed_h_and_matching(self):
        payload=propose(dress_polynomial(quadratic_fixture(),((0,1),),4),4,2)
        bad=copy.deepcopy(payload)
        bad['certificate']['hamiltonian'][0]['coefficient']='4'
        with self.assertRaises(ValueError): replay(bad)
        for matching in ([[0,2]],[[0,1],[1,2]],[[True,1]]):
            bad=copy.deepcopy(payload);bad['matching']=matching
            with self.assertRaises(ValueError): replay(bad)

    def test_refuses_bad_occupied_witness(self):
        payload=propose(dress_polynomial(quadratic_fixture(),((0,1),),4),4,2)
        for columns in ([["1","0"],["2","0"],["3","0"],["4","0"]],[[1,0]],[[1.0,0]]*4):
            bad=copy.deepcopy(payload);bad['occupied_columns']=columns
            with self.assertRaises(ValueError): replay(bad)

    def test_rounding_and_factor_corruption_paid_by_exact_residual(self):
        h=dress_polynomial(quadratic_fixture(),((0,1),),4)
        payload=propose(h,4,2,denominator=10**4)
        result=replay(payload)
        self.assertGreater(F(result['lower_replay']['residual_l1']),0)
        bad=copy.deepcopy(payload)
        bad['certificate']['blocks'][0]['factor'][0][0]+=10000
        worse=replay(bad)
        self.assertLessEqual(worse['lower_float'], result['upper_float'])
        self.assertGreater(F(worse['lower_replay']['residual_l1']),F(result['lower_replay']['residual_l1']))

if __name__=='__main__': unittest.main()
