import unittest
from copy import deepcopy
from fractions import Fraction as F
from experiments.v7_certificate import Generator, Piece, derive_certificate, check_certificate, residual_records, norm_witness

class CertificateTests(unittest.TestCase):
    def test_rotation_residual_and_dense_error(self):
        # H=Z/2, G(X)=-Y, G(Y)=X. First-order P=X-tY, residual=-tX.
        g=lambda:Generator({'Z':F(1,2)},F(0),1)
        ps=[Piece(F(1,10),({'X':F(1)},{'Y':F(-1)}))]
        records=residual_records(g(), {'X':F(1)}, ps)
        self.assertEqual(records[-1][1], {'X':F(1)})
        w=derive_certificate(g(), {'X':F(1)}, ps)
        ans=check_certificate(g(), {'X':F(1)}, ps,w,F(1,200),expected_time=F(1,10))
        self.assertEqual(ans['status'],'certified')
        self.assertEqual(F(ans['bound']),F(1,200))
        import math
        actual=math.hypot(math.cos(.1)-1, -math.sin(.1)+.1)
        self.assertLessEqual(actual,float(F(ans['bound'])))

    def test_initial_error_and_jump_are_charged(self):
        g=lambda:Generator({},F(0),1)
        ps=[Piece(F(1),({'X':F(2)},)),Piece(F(1),({'X':F(5)},))]
        w=derive_certificate(g(),{'X':F(1)},ps)
        a=check_certificate(g(),{'X':F(1)},ps,w,F(3),expected_time=F(2))
        self.assertEqual(a['status'],'over_tolerance')
        self.assertEqual(F(a['bound']),4)

    def test_damping_and_anticommuting_identity(self):
        g=Generator({},F(1,5),2)
        self.assertEqual(g.apply({'XY':F(3)}),{'XY':F(-6,5)})
        w,c=norm_witness({'X':F(3),'Y':F(4)},'firstfit')
        self.assertEqual(w,[{'labels':['X','Y'],'upper':'5'}])
        self.assertEqual(c['group_comparisons'],1)

    def test_tampering_fails(self):
        g=lambda:Generator({'Z':F(1,2)},F(0),1)
        ps=[Piece(F(1),({'X':F(1),'Y':F(1)},))]
        w=derive_certificate(g(),{},ps,'weighted')
        for field in ('coverage','bound','sqrt','duplicate'):
            bad=deepcopy(w)
            if field=='coverage': del bad['witnesses']['jump:0']
            elif field=='bound': bad['claimed_bound']='0'
            elif field=='sqrt': bad['witnesses']['jump:0'][0]['upper']='0'
            else: bad['witnesses']['jump:0'][0]['labels'].append('X')
            self.assertEqual(check_certificate(g(),{},ps,bad,F(10),expected_time=F(1))['status'],'rejected')
        # A commuting pair has no l2 norm certificate.
        g2=lambda:Generator({},F(0),2)
        ps2=[Piece(F(1),({'XI':F(1),'IX':F(1)},))]
        w2=derive_certificate(g2(),{},ps2)
        w2['witnesses']['jump:0']=[dict(labels=['XI','IX'],upper='2')]
        self.assertEqual(check_certificate(g2(),{},ps2,w2,F(10),expected_time=F(1))['status'],'rejected')

    def test_horizon_and_parser(self):
        g=lambda:Generator({},F(0),1)
        ps=[Piece(F(1),({'X':F(1)},))]
        w=derive_certificate(g(),{},ps)
        self.assertEqual(check_certificate(g(),{},ps,w,F(2),expected_time=F(2))['status'],'rejected')
        for value in (1.0, '9'*6000, '1/0', '1e10'):
            bad=deepcopy(w); bad['witnesses']['jump:0'][0]['upper']=value
            self.assertEqual(check_certificate(g(),{},ps,bad,F(2),expected_time=F(1))['status'],'rejected')

    def test_bernstein_cancellation_and_checker(self):
        # G=0, P(t)=(t-t^2/2)X, P'= (1-t)X: exact integral 1/2.
        g=lambda:Generator({},F(0),1)
        ps=[Piece(F(1),({}, {'X':F(1)}, {'X':F(-1,2)}))]
        power=derive_certificate(g(),{},ps,integration_basis='power')
        bern=derive_certificate(g(),{},ps,integration_basis='bernstein')
        self.assertEqual(F(power['claimed_bound']),F(3,2))
        self.assertEqual(F(bern['claimed_bound']),F(1,2))
        self.assertEqual(check_certificate(g(),{},ps,bern,F(1,2),expected_time=F(1))['status'],'certified')
        bad=deepcopy(bern); bad['integration_basis']='sample_endpoints'
        self.assertEqual(check_certificate(g(),{},ps,bad,F(1),expected_time=F(1))['status'],'rejected')

    def test_invalid_model_and_budget(self):
        for h,gamma,n in [({'X':.5},F(0),1),({},F(-1),1),({},F(0),33)]:
            with self.assertRaises(ValueError): Generator(h,gamma,n)
        g=lambda:Generator({},F(0),1)
        ps=[Piece(F(-1),({'X':F(1)},))]
        self.assertEqual(check_certificate(g(),{},ps,{},F(1),expected_time=F(1))['status'],'rejected')

if __name__=='__main__': unittest.main()
