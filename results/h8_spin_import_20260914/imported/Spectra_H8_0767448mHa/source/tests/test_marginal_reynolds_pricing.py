import contextlib,io,unittest
from fractions import Fraction as F
import numpy as np
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_moment_pricing import MomentDictionary
from experiments.marginal_reynolds_pricing import ReynoldsMomentDictionary,complement_polynomial


class ReynoldsPricingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw=build(4,4,F(1));cls.data={k:raw[k] for k in ('modes','particles','hamiltonian')}
        with contextlib.redirect_stdout(io.StringIO()):
            cls.base=MomentDictionary(cls.data,'-14',feature_degree=2)
            cls.reduced=ReynoldsMomentDictionary(cls.data,'-14')

    def test_full_group_columns_metric_and_dual_pricing(self):
        b,d=self.base,self.reduced;L=d.reynolds_L
        self.assertLess(d.qrows,b.qrows)
        self.assertEqual(d.metric.shape,(2*d.qrows+1,len(d.orbits)))
        for orbit in d.orbits:self.assertTrue(all(sum(p)%2==0 for p in orbit))
        expected=b.metric[:,d.retained_feature_indices].toarray()
        np.testing.assert_allclose(d.metric.toarray(),np.vstack((L@expected[:b.qrows],L@expected[b.qrows:2*b.qrows],expected[-1:])),atol=1e-10)
        for label in d.labels(3)[::7]:
            members=d.atom_members(label)
            for f,r,o in members:self.assertEqual(o&~r,0)
            projected=[L@b.column(member) for member in members]
            for column in projected:np.testing.assert_allclose(d.column(label),column,atol=1e-10)
        dual=np.random.default_rng(7).normal(size=2*d.qrows+1)
        selected,maximum,_=d.price(dual,set(),4)
        self.assertGreater(maximum,0)
        for block,label in selected:self.assertGreater(-dual[block*d.qrows:(block+1)*d.qrows]@d.column(label),0)
        self.assertEqual(d.qrows,len(L))

    def test_exact_refusals_and_particle_hole_identity(self):
        d=self.reduced
        self.assertEqual(complement_polynomial(complement_polynomial({3:F(2),4:F(-1)})),{3:F(2),4:F(-1)})
        with self.assertRaisesRegex(ValueError,'Reflection'):d.check_invariant({1:1,2:1})
        with self.assertRaises(ValueError):d.atom_members(('positive',1,2))
        with self.assertRaises(ValueError):ReynoldsMomentDictionary(dict(self.data,particles=2),'-14')
        with self.assertRaises(ValueError):ReynoldsMomentDictionary(self.data,'-14',reflection=1)


if __name__=='__main__':unittest.main()
