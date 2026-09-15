import copy
import unittest
from fractions import Fraction as F
from math import sqrt
from experiments.v8_fractional_cover import propose_fractional_cover, check_fractional_cover


class FractionalCoverTests(unittest.TestCase):
    def test_odd_cycle_strictly_beats_every_partition(self):
        # Anticommutation graph C5: any disjoint partition uses >=1 singleton
        # and at most 2 edges, hence optimum 1+2sqrt(2). Five half-edges beat it.
        op={p:F(1) for p in ('IX','IY','XI','YX','ZY')}
        edges=[['IX','IY'],['IX','ZY'],['IY','YX'],['XI','YX'],['XI','ZY']]
        w,_=propose_fractional_cover(op,groups=edges,iterations=4)
        self.assertEqual(check_fractional_cover(op,w)['status'],'certified')
        self.assertAlmostEqual(float(F(w['claimed_bound'])),5/sqrt(2))
        self.assertLess(float(F(w['claimed_bound'])),1+2*sqrt(2))
        generated,_=propose_fractional_cover(op,iterations=4)
        self.assertLess(float(F(generated['claimed_bound'])),1+2*sqrt(2))

    def test_full_clique_partition_is_fairly_available(self):
        op={'X':F(1),'Y':F(1),'Z':F(1)}
        w,_=propose_fractional_cover(op,groups=[['X','Y'],['Y','Z'],['X','Z']])
        self.assertAlmostEqual(float(F(w['claimed_bound'])),sqrt(3))

    def test_signed_rational_reconstruction(self):
        op={'IX':F(-7,3),'IY':F(2,5),'XI':F(1,9),'YX':F(-3),'ZY':F(4)}
        w,_=propose_fractional_cover(op,iterations=6)
        self.assertEqual(check_fractional_cover(op,w)['status'],'certified')
        w['groups'][0]['coefficients'][next(iter(w['groups'][0]['coefficients']))]='99'
        self.assertEqual(check_fractional_cover(op,w)['status'],'rejected')

    def test_rounding_and_extra_coefficients_rejected(self):
        op={'X':F(1),'Y':F(1)}
        w,_=propose_fractional_cover(op)
        for upper in ('-2','7/5',1.42,'1'*5003):
            bad=copy.deepcopy(w);bad['groups'][0]['upper']=upper
            self.assertEqual(check_fractional_cover(op,bad)['status'],'rejected')
        bad=copy.deepcopy(w);bad['groups'][0]['coefficients']['I']='1'
        self.assertEqual(check_fractional_cover(op,bad)['status'],'rejected')

    def test_nonanti_and_missing_fields_fail_closed(self):
        bad=dict(schema='v8-cover-2',groups=[dict(coefficients={'IX':'1','XI':'1'},upper='2')],claimed_bound='2')
        self.assertEqual(check_fractional_cover({'IX':1,'XI':1},bad)['status'],'rejected')
        for bad in ({},None,{'schema':'v8-cover-2','groups':[None]}):
            self.assertEqual(check_fractional_cover({'X':1},bad)['status'],'rejected')

    def test_budgets_and_labels(self):
        for op,kwargs in (({'X':1,'XX':2},{}),({'Z':1.0},{}),({'X':1},{'iterations':129}),({'X':1},{'groups':[['X']]*129}),({'X':1},{'groups':[['X','X']]})):
            with self.assertRaises(ValueError):propose_fractional_cover(op,**kwargs)
        w=dict(schema='v8-cover-2',groups=[dict(coefficients={'X':'1'},upper='1')]*513,claimed_bound='513')
        self.assertEqual(check_fractional_cover({'X':513},w)['status'],'rejected')

    def test_empty_zero(self):
        w,_=propose_fractional_cover({})
        self.assertEqual(check_fractional_cover({},w)['bound'],'0')

    def test_full_evolution_and_tampering(self):
        from experiments.v8_fractional_cover import adaptive_cover, check_evolution_cover
        from experiments.v7_headroom import model,TOL
        from experiments.v7_certificate import Generator
        h,o=model(3,'xxz')
        for enabled,order in ((False,12),(True,11)):
            p,w,_=adaptive_cover(Generator(h,F(2),3,512),o,F(1,2),TOL,enabled=enabled)
            self.assertEqual(len(p.coefficients)-1,order)
            def check(witness,horizon=F(1,2)):
                return check_evolution_cover(Generator(h,F(2),3,512),o,[p],witness,TOL,expected_time=horizon)
            self.assertEqual(check(w)['status'],'certified')
            self.assertEqual(check(w,F(1,3))['status'],'rejected')
            bad=copy.deepcopy(w);bad['witnesses']['residual:0:0']['claimed_bound']='1'
            self.assertEqual(check(bad)['status'],'rejected')
            bad=copy.deepcopy(w);bad['claimed_bound']='0'
            self.assertEqual(check(bad)['status'],'rejected')

    def test_large_partition_budget_matches_v7(self):
        from experiments.v8_fractional_cover import adaptive_cover, check_evolution_cover
        from experiments.v7_headroom import model,TOL
        from experiments.v7_certificate import Generator
        h,o=model(6,'xxz')
        p,w,_=adaptive_cover(Generator(h,F(0),6,512),o,F(1,5),TOL)
        result=check_evolution_cover(Generator(h,F(0),6,512),o,[p],w,TOL,expected_time=F(1,5))
        self.assertEqual(result['status'],'certified')

if __name__=='__main__': unittest.main()
