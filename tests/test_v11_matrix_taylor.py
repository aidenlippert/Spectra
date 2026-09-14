import unittest
from fractions import Fraction as F
from itertools import product
import numpy as np
from experiments.v11_matrix_taylor import MatrixAction,matrix_taylor
from experiments.v7_certificate import Generator,check_certificate,BudgetExceeded
from experiments.v8_integer_taylor import fraction_free_taylor
from experiments.v7_headroom import model,TOL


class MatrixTaylorTests(unittest.TestCase):
    def test_every_two_qubit_pauli_column_matches_original(self):
        h={'ZI':F(1,2),'IX':F(2,3),'XY':F(-3,5),'YZ':F(2,7),'II':F(9)}
        gamma=F(2,5);engine=MatrixAction(h,gamma,2);gen=Generator(h,gamma,2)
        for p in map(''.join,product('IXYZ',repeat=2)):
            pair,q0=engine.initial({p:F(1)})
            self.assertEqual(q0,1)
            self.assertEqual(engine.to_pauli(pair),{p:1})
            got=engine.to_pauli(engine.apply(pair))
            self.assertEqual(got,{q:int(c*engine.q) for q,c in gen.apply({p:F(1)}).items()})

    def test_dense_physical_action_independent(self):
        mats={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),
              'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1,-1])}
        dense=lambda p:np.kron(mats[p[0]],mats[p[1]])
        h={'ZI':F(1),'YX':F(-2),'XY':F(3)};gamma=F(4)
        o={'XZ':F(5),'YI':F(-7),'II':F(3)}
        e=MatrixAction(h,gamma,2);pair,q0=e.initial(o)
        A=sum(float(c)*dense(p) for p,c in o.items())
        H=sum(float(c)*dense(p) for p,c in h.items())
        expected=1j*(H@A-A@H)
        for p in ['XI','YI','ZI','IX','IY','IZ']:
            P=dense(p);expected+=P@A@P-A
        r,im=e.apply(pair)
        np.testing.assert_array_equal(r+1j*im,expected*e.q*q0)

    def test_exact_polynomial_and_certificate_equality(self):
        for n,gamma,T in [(3,F(2),F(1,2)),(3,F(1,5),F(1,5)),(4,F(0),F(1,5))]:
            h,o=model(n,'xxz')
            p,w,c=matrix_taylor(Generator(h,gamma,n,512),o,T,TOL)
            original,ow,_=fraction_free_taylor(Generator(h,gamma,n,512),o,T,TOL)
            self.assertEqual(p,original);self.assertEqual(w['witnesses'],ow['witnesses'])
            self.assertEqual(w['claimed_bound'],ow['claimed_bound'])
            result=check_certificate(Generator(h,gamma,n,512),o,[p],w,TOL,expected_time=T)
            self.assertEqual(result['status'],'certified')

    def test_rational_initial_and_stationary(self):
        p,w,_=matrix_taylor(Generator({},0,1),{'X':F(1,3)},F(1),F(0))
        self.assertEqual(p.coefficients,({'X':F(1,3)},));self.assertEqual(w['claimed_bound'],'0')
        e=MatrixAction({'Z':F(1)},F(0),1)
        pair,q0=e.initial({'Y':F(2,3),'X':F(1,5)})
        self.assertEqual(q0,15);self.assertEqual(e.to_pauli(pair),{'X':3,'Y':10})

    def test_pre_kernel_overflow_falls_back_exactly(self):
        e=MatrixAction({'Z':F(1)},F(0),1)
        pair,_=e.initial({'X':F(1<<62)})
        self.assertEqual(pair[0].dtype,np.dtype('int64'))
        result=e.apply(pair)
        self.assertEqual(result[0].dtype,np.dtype('O'))
        self.assertEqual(e.to_pauli(result),{'Y':-(1<<63)})
        self.assertGreater(e.cost['object_kernels'],0)
        pair,_=e.initial({'Y':F(1<<90)})
        self.assertEqual(e.to_pauli(e.apply(pair)),{'X':1<<91})

    def test_refusals_and_no_truncation(self):
        with self.assertRaises(ValueError):MatrixAction({},F(-1),1)
        with self.assertRaises(BudgetExceeded):MatrixAction({},0,7)
        with self.assertRaises(ValueError):MatrixAction({'Q':F(1)},0,1)
        e=MatrixAction({'X':F(1),'Z':F(1)},0,1,cap=1)
        pair,_=e.initial({'Y':F(1)})
        with self.assertRaises(BudgetExceeded):e.to_pauli(e.apply(pair))
        e=MatrixAction({},0,1)
        with self.assertRaises(ValueError):e.to_pauli((np.diag([1,0]),np.zeros((2,2),dtype=np.int64)))
        h,o=model(3,'xxz')
        with self.assertRaises(BudgetExceeded):matrix_taylor(Generator(h,2,3,512),o,F(1,2),TOL,max_order=0)


if __name__=='__main__':unittest.main()
