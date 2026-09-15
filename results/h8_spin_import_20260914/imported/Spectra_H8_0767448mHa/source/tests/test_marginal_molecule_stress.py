import json, unittest
from pathlib import Path
from fractions import Fraction as F
from experiments.marginal_symbolic import decode
from experiments.marginal_transfer import upper_ed
from experiments.marginal_transfer_verify import replay
from experiments.marginal_molecule_stress import parity_generators

class StressFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=json.loads((Path(__file__).resolve().parents[1]/'results/marginal_molecule_stress/h4_square_sto3g.json').read_text())
        cls.h=decode(cls.d['hamiltonian'],8,4)
    def test_common_denominator_and_exact_l1_budget(self):
        self.assertEqual(self.d['coefficient_denominator'],10**12)
        self.assertTrue(all(10**12 % c.denominator==0 for c in self.h.values()))
        err=F(self.d['numerical_integral_perturbation_l1_upper'])
        self.assertGreaterEqual(err,0); self.assertLess(err,F(1,10**8))
    def test_integer_upper_witness_replays(self):
        w=upper_ed(self.h,8,4)
        self.assertTrue(w['amplitudes'] and all(type(x) is int for x in w['amplitudes']))
        c={'modes':8,'particles':4,'hamiltonian':self.d['hamiltonian'],
           'number_multiplier':[],'b':'-10000','blocks':[], 'denominator':1,'independent_upper':w}
        independently=replay(c)
        self.assertEqual(F(w['upper']),F(independently['upper']))
        w['upper']='1000'
        self.assertEqual(F(independently['upper']),F(replay(c)['upper']))

    def test_binary_nullspace_matches_exhaustive_symmetry_scan(self):
        generators=parity_generators(self.h,8);span={0}
        for g in generators:span|={v^g for v in list(span)}
        expected={s for s in range(256) if all(sum((s>>i)&1 for c,i in w)%2==0 for w in self.h)}
        self.assertEqual(span,expected)
        self.assertEqual(len(span),2**len(generators))
        self.assertEqual(len(generators),3)

if __name__=='__main__': unittest.main()
