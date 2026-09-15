import unittest
from fractions import Fraction as F
from experiments.v7_certificate import Generator,check_certificate
from experiments.v9_polynomial import krylov_basis,collocation_piece,projected_taylor_piece,best_witness

class PolynomialTests(unittest.TestCase):
    def test_rotation_and_exact_initial(self):
        h={'Z':F(1,2)};o={'X':F(1)};T=F(1,2);eps=F(1,1000)
        basis=krylov_basis(Generator(h,F(0),1,512),o,4)
        self.assertEqual(basis['A'].shape,(2,2))
        for construct in (collocation_piece,projected_taylor_piece):
            piece,cost=construct(basis,o,T,6)
            self.assertEqual(piece.coefficients[0],o)
            w,_=best_witness(Generator(h,F(0),1,512),o,piece,eps)
            r=check_certificate(Generator(h,F(0),1,512),o,[piece],w,eps,expected_time=T)
            self.assertEqual(r['status'],'certified',r)
            bad=check_certificate(Generator(h,F(0),1,512),o,[piece],w,eps,expected_time=F(1))
            self.assertEqual(bad['status'],'rejected')

    def test_two_qubit_interacting(self):
        h={'ZI':F(1,3),'XX':F(2,5),'IZ':F(1,7)};o={'ZI':F(1)};T=F(1,5);eps=F(1,1000)
        basis=krylov_basis(Generator(h,F(1,5),2,512),o,8)
        piece,cost=collocation_piece(basis,o,T,6)
        w,_=best_witness(Generator(h,F(1,5),2,512),o,piece,eps)
        r=check_certificate(Generator(h,F(1,5),2,512),o,[piece],w,eps,expected_time=T)
        self.assertEqual(r['status'],'certified',r)
        self.assertGreater(cost['linear_system_dimension'],0)

    def test_request_caps(self):
        g=Generator({'Z':F(1)},F(0),1,512)
        with self.assertRaises(ValueError):krylov_basis(g,{'X':F(1)},33)
        basis=krylov_basis(g,{'X':F(1)},2)
        with self.assertRaises(ValueError):collocation_piece(basis,{'X':F(1)},F(1),13)
        with self.assertRaises(ValueError):collocation_piece(basis,{'X':F(1)},F(0),4)

if __name__=='__main__':unittest.main()
