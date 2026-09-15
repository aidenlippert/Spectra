import contextlib,io,unittest
from fractions import Fraction as F
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_moment_pricing import MomentDictionary
from experiments.marginal_reynolds_pricing import ReynoldsMomentDictionary
from experiments.marginal_conditional_coordinates import conditional_map


class ConditionalCoordinateTests(unittest.TestCase):
    def test_binomial_global_means_and_invariant_rank(self):
        raw=build(8,4,F(1));data={k:raw[k] for k in ('modes','particles','hamiltonian')}
        with contextlib.redirect_stdout(io.StringIO()):
            base=MomentDictionary(data,'-19',feature_degree=0,proof_degree=4)
            T,stats=conditional_map(base)
            self.assertAlmostEqual(float((T@base.column(('positive',1,1)))[0]),.5)
            self.assertAlmostEqual(float((T@base.column(('charge',0,0)))[0]),float(F(71,69)))
            self.assertLess(stats['maximum_direct_count_error'],1e-9)
            a=ReynoldsMomentDictionary(data,'-19',feature_degree=0,proof_degree=4)
            b=ReynoldsMomentDictionary(data,'-19',feature_degree=0,proof_degree=4,coordinates='conditional')
            self.assertEqual(a.qrows,b.qrows)
            self.assertEqual(a.metric.shape,b.metric.shape)
            self.assertEqual(b.stats['reynolds_coordinates'],'conditional')

    def test_rank_loss_and_unknown_coordinates_refused(self):
        raw=build(4,4,F(1));data={k:raw[k] for k in ('modes','particles','hamiltonian')}
        with self.assertRaisesRegex(ValueError,'coordinates'):ReynoldsMomentDictionary(data,'-14',coordinates='other')
        with contextlib.redirect_stdout(io.StringIO()),self.assertRaises(ValueError):
            ReynoldsMomentDictionary(data,'-14',coordinates='conditional')


if __name__=='__main__':unittest.main()
