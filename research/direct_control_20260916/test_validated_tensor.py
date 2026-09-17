import unittest
from fractions import Fraction as F
import numpy as np
from .tensor_operator import build,element
from .validated_tensor import zip_apply,load_mps,expectation,opnorm,identity_mpo


def dense(a):
    x=a[0][0]
    for t in a[1:]:x=np.tensordot(x,t,axes=([-1],[0]))
    return x.reshape(-1)


class TestValidatedTensor(unittest.TestCase):
    def test_compression_and_action_enclose_independent_dense_oracle(self):
        rng=np.random.default_rng(710)
        ts=[rng.normal(size=(1,2,2))+1j*rng.normal(size=(1,2,2)),
            rng.normal(size=(2,2,3))+1j*rng.normal(size=(2,2,3)),
            rng.normal(size=(3,2,2))+1j*rng.normal(size=(3,2,2)),
            rng.normal(size=(2,2,1))+1j*rng.normal(size=(2,2,1))]
        d={'modes':4,'hamiltonian':[{'word':[], 'coefficient':'3/7'},
              {'word':[(1,0),(0,3)],'coefficient':'2/3'},
              {'word':[(1,3),(0,0)],'coefficient':'2/3'},
              {'word':[(1,2),(0,2)],'coefficient':'-4/5'}]}
        op=build(d)
        labels=[int(format(x,'04b')[::-1],2) for x in range(16)]
        h=np.array([[float(element(op,i,j)) for j in labels] for i in labels])
        for bond in (1,2,4):
            out,e,_=zip_apply(ts,op,bond)
            actual=np.linalg.norm(h@dense(ts)-dense(out))
            self.assertLessEqual(actual,e)
            replay,er,_=zip_apply(ts,op,bond,proposal=out)
            self.assertGreaterEqual(er,actual)
            self.assertTrue(all(np.array_equal(x,y) for x,y in zip(out,replay)))
        out,e,_=zip_apply(ts,op,16)
        self.assertLess(e,1e-6)

    def test_expectation_enclosure(self):
        rng=np.random.default_rng(99)
        ts=[rng.normal(size=(1,2,2))+1j*rng.normal(size=(1,2,2)),
            rng.normal(size=(2,2,1))+1j*rng.normal(size=(2,2,1))]
        v=dense(ts)
        for diag,exact in [(None,np.vdot(v,v)),({0:(1,-1)},np.vdot(v,np.diag([1,1,-1,-1])@v))]:
            c,e=expectation(ts,diag)
            self.assertLessEqual(abs(c-exact),e)

    def test_corrupt_proposal_charged_and_bad_arithmetic_refused(self):
        ts=[np.array([[[1.],[2.]]],complex),np.array([[[3.],[4.]]],complex)]
        op=identity_mpo(2);out,e,_=zip_apply(ts,op,2)
        bad=[x.copy() for x in out];bad[-1][0,0,0]+=1
        _,eb,_=zip_apply(ts,op,2,proposal=bad)
        self.assertGreater(eb,.1)
        with self.assertRaises(ArithmeticError):opnorm(np.array([[float('nan')]]))


if __name__=='__main__':unittest.main()
