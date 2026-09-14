import copy
from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import mono,add,scale,product,encode
from research.collective_completion_20260914.spin_replay import check,setup,alpha_shift

class SpinReplayTests(unittest.TestCase):
    def fixture(self):
        h=add(*(mono(((1,i),(0,i)),F(1,3)) for i in range(4)))
        data={'modes':4,'particles':2,'hamiltonian':encode(h)}
        core={**data,'operator_degree':3,'number_multiplier':encode(mono((),F(1,3))),
              'b':'2/3','denominator':1,'blocks':[]}
        cert={**data,'kind':'balanced_spin_sos_v1','alpha_multiplier':[],'core':core}
        return data,cert
    def test_exact_endpoint(self):
        d,c=self.fixture();self.assertEqual(F(check(d,c)['lower']),F(2,3))
    def test_spin_defect_charged(self):
        d,c=self.fixture();h=setup(d)[0]
        z=add(mono(((1,0),(0,0))),mono(((1,1),(0,1)),-1))
        d['hamiltonian']=encode(add(h,scale(z,F(1,1000))));c['hamiltonian']=d['hamiltonian']
        self.assertEqual(F(check(d,c)['lower']),F(2,3)-F(1,500))
    def test_false_balanced_only_bound_refused(self):
        d,c=self.fixture();h=scale(product(alpha_shift(4,2),alpha_shift(4,2)),-1)
        d['hamiltonian']=encode(h);c['hamiltonian']=d['hamiltonian']
        c['alpha_multiplier']=encode(scale(alpha_shift(4,2),-1));c['core']['hamiltonian']=[]
        c['core']['b']='0';c['core']['number_multiplier']=[]
        with self.assertRaisesRegex(ValueError,'Auxiliary'):check(d,c)
    def test_corruption_refusals(self):
        d,c=self.fixture()
        for key,value in [('particles',1),('modes',3),('kind','unknown')]:
            bad=copy.deepcopy(c);bad[key]=value
            with self.assertRaises(ValueError):check(d,bad)
        bad=copy.deepcopy(c);bad['alpha_multiplier']=encode(mono(((1,0),(0,1))))
        with self.assertRaises(ValueError):check(d,bad)
        bad=copy.deepcopy(c);bad['core']['hamiltonian']=[]
        with self.assertRaises(ValueError):check(d,bad)
if __name__=='__main__':unittest.main()
