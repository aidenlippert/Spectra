import copy
from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import mono,add,encode,scale,product,canonical,adj
from research.collective_completion_20260914.spin_screen import spin_squared,ladder_ideal,check,check_sector
from research.collective_completion_20260914.spin_replay import alpha_shift

class SpinScreenTests(unittest.TestCase):
    def pair(self):
        data={'modes':4,'particles':2,'hamiltonian':encode(add(*(mono(((1,i),(0,i)),F(1,3)) for i in range(4))))}
        core={**data,'operator_degree':3,'number_multiplier':encode(mono((),F(1,3))),'b':'2/3','denominator':1,'blocks':[]}
        c={**data,'kind':'spin_sector_sos_v1','alpha_multiplier':[],'core':core,'casimir_multiplier':'0','magnetization':0,'singlet':True}
        t=copy.deepcopy(c);t.update(magnetization=1,singlet=False)
        return data,c,t
    def test_complete_constant_energy(self):
        d,s,t=self.pair();self.assertEqual(F(check(d,s,t)['lower']),F(2,3))
    def test_two_sector_minimum_and_refusals(self):
        d,s,t=self.pair();t['core']['b']='1/2';self.assertLessEqual(F(check(d,s,t)['lower']),F(1,2))
        with self.assertRaises(ValueError):check(d,s,s)
        t['casimir_multiplier']='1'
        with self.assertRaises(ValueError):check(d,s,t)
    def test_casimir_singlet_identity(self):
        d,s,t=self.pair();d['hamiltonian']=encode(spin_squared(4));s['hamiltonian']=d['hamiltonian']
        s['casimir_multiplier']='1';s['core']['hamiltonian']=[];s['core']['number_multiplier']=[];s['core']['b']='0'
        self.assertEqual(F(check_sector(d,s)['lower']),0)
        s['singlet']=False
        with self.assertRaises(ValueError):check_sector(d,s)
    def test_spin_ladder_zero_compression_independently(self):
        from research.collective_completion_20260914.angles import act
        singlets=[{3:F(1)},{12:F(1)},{9:F(1),6:F(-1)}]
        def matrix_element(left,p,right):
            out=F(0)
            for x,a in right.items():
                for w,c in p.items():
                    y,s=act(w,x)
                    if s:out+=left.get(y,0)*s*c*a
            return out
        for i in (1,3):
            for j in (0,2):
                p=ladder_ideal(4,mono(((1,i),(0,j))))
                for left in singlets:
                    for right in singlets:self.assertEqual(matrix_element(left,p,right),0)
        for v in singlets:self.assertEqual(matrix_element(v,spin_squared(4),v),0)
        self.assertEqual(matrix_element({5:F(1)},spin_squared(4),{5:F(1)}),2)
    def test_ladder_mutation_refused(self):
        d,s,t=self.pair();s['spin_ladder_multiplier']=encode(mono(((1,0),(0,1))))
        with self.assertRaises(ValueError):check_sector(d,s)
        s['spin_ladder_multiplier']=encode(mono(((1,1),(0,0))))
        with self.assertRaisesRegex(ValueError,'Auxiliary'):check_sector(d,s)
if __name__=='__main__':unittest.main()
