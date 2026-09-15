import copy
from fractions import Fraction as F
import math
import random
import unittest
from research.correlated_pair_20260913.test_mps import correlated_fixture
from research.correlated_pair_20260913.mps_exact import check as exact_check
from research.reconstruction_compression_20260914.upper_interval import plus,times,rational_interval,check

class IntervalTests(unittest.TestCase):
    def assertContains(self,interval,value):
        self.assertLessEqual(F(float(interval[0])),value)
        self.assertGreaterEqual(F(float(interval[1])),value)

    def test_primitive_enclosures_against_rationals(self):
        rng=random.Random(432)
        special=[0.,1.,-1.,math.ldexp(1.,-1022),math.ldexp(1.,-1074),1e100,-1e100]
        vals=special+[math.ldexp(rng.uniform(-1,1),rng.randint(-400,400)) for _ in range(200)]
        for a,b in zip(vals,vals[3:]+vals[:3]):
            self.assertContains(plus(a,a,b,b),F(a)+F(b))
            self.assertContains(times(a,a,b,b),F(a)*F(b))
        for _ in range(100):
            q=F(rng.randint(-10**30,10**30),rng.randint(1,10**30))
            self.assertContains(rational_interval(q),q)

    def test_exact_energy_and_strict_refusals(self):
        data,state=correlated_fixture();exact=exact_check(data,state);u=F(exact['upper_Ha'])
        rec=check(data,state,u+F(1,10**10));self.assertEqual(rec['status'],'certified_upper')
        norm=list(map(F,rec['norm_enclosure']));self.assertLessEqual(norm[0],5);self.assertGreaterEqual(norm[1],5)
        self.assertEqual(check(data,state,u-F(1,10**10))['status'],'endpoint_refuted')
        self.assertEqual(check(data,state,u)['status'],'ambiguous_refused')
        self.assertEqual(check(data,state,u,True)['status'],'certified_upper_exact_fallback')

    def test_invalid_witnesses(self):
        data,base=correlated_fixture()
        for change in ('charge','boundary','binding','zero','duplicate','integer','denominator'):
            state=copy.deepcopy(base)
            if change=='charge':state['tensors'][0][0][1]=1
            if change=='boundary':state['bond_charges'][-1]=[[2,0]]
            if change=='binding':state['fixture_sha256']='wrong'
            if change=='zero':state['tensors'][1]=[]
            if change=='duplicate':state['tensors'][0].append(state['tensors'][0][0])
            if change=='integer':state['tensors'][0][0][-1]=0.5
            if change=='denominator':state['denominator']=-1
            with self.subTest(change=change),self.assertRaises(ValueError):check(data,state,F(100))

if __name__=='__main__':unittest.main()
