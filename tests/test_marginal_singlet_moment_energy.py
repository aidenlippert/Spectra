import copy,json,unittest
from pathlib import Path
from results.marginal_graded_hubbard8.discovery.singlet_moment_energy import replay,polynomial_upper,upper_monomials
from experiments.marginal_symmetry_moments import projected_moments
ROOT=Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8'
class SingletMomentEnergyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.c=json.loads((ROOT/'singlet_moment_energy/certificate.json').read_text())
    def test_exact_moment_energy_interval(self):
        r=replay(self.c)
        self.assertEqual(r['lower'],'-106/25');self.assertLess(r['upper_float'],-4.2356)
        self.assertLess(r['width_float'],.0044)
        self.assertLess(r['moment_work']['unique_determinant_source_actions'],1200)
    def test_refuse_unproved_lower_and_bad_upper(self):
        c=copy.deepcopy(self.c);c['lower']='-4.23'
        with self.assertRaisesRegex(ValueError,'positive moment response'):replay(c)
        c=copy.deepcopy(self.c);c['upper_polynomial_coefficients']=[[0]*14]
        with self.assertRaisesRegex(ValueError,'nonzero upper'):replay(c)
        c=copy.deepcopy(self.c);c['upper_polynomial_coefficients'][0][0]=10**16
        with self.assertRaisesRegex(ValueError,'bounded integer upper'):replay(c)
    def test_recurrence_and_upper_exact_scalar_sanity(self):
        # Nonorthogonal scalar renewal with G=2 and moments3,7,20.
        from fractions import Fraction as F
        self.assertEqual(projected_moments([[[2]],[[3]],[[7]],[[20]]]),[[[F(5,2)]],[[F(23,4)]]])
        # Direct polynomial Rayleigh check for H=diag(1,3), V=(1,1).
        upper,norm=polynomial_upper([[[2]],[[4]],[[10]],[[28]]],[[1],[1]])
        self.assertEqual(norm,20);self.assertEqual(upper,F(13,5))
        with self.assertRaises(ValueError):projected_moments([[[0]],[[1]],[[2]]])
        coefficients=upper_monomials([[1],[2],[3]])
        moments=[[[1+3**k]] for k in range(6)]
        upper,norm=polynomial_upper(moments,coefficients)
        def value(x):
            z=F(x-8,12);return 1+2*z+3*(2*z*z-1)
        expected_norm=value(1)**2+value(3)**2
        self.assertEqual(norm,expected_norm)
        self.assertEqual(upper,(value(1)**2+3*value(3)**2)/expected_norm)
    def test_second_response_block_and_fixed_moment_budget(self):
        c=copy.deepcopy(self.c);p=c.pop('response_polynomial');c['response_polynomials']=[p,['1']]
        r=replay(c);self.assertEqual(r['response_dimension'],28)
        c['response_polynomials']=[p,p]
        with self.assertRaises(ValueError):replay(c)
        c=copy.deepcopy(self.c);c['response_polynomial']=['1']*12
        with self.assertRaisesRegex(ValueError,'moment-order24'):replay(c)
        c=copy.deepcopy(self.c);c['upper_polynomial_coefficients']=[[1]*14]*13
        with self.assertRaisesRegex(ValueError,'moment-order24'):replay(c)
if __name__=='__main__':unittest.main()
