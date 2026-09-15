import copy
from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import mono,add,adj,canonical,encode
from research.correlated_pair_20260913.density_conditioned.exact_density import allocation,closure_demo
from research.correlated_pair_20260913.self_consistent.fixed_guide import build,check
from research.correlated_pair_20260913.complement import ldlt_positive


class AlgebraChecks(unittest.TestCase):
    def test_flow_and_explicit_cut(self):
        self.assertTrue(allocation([F(1),F(2)],{(0,1):F(-2)})['feasible'])
        bad=allocation([F(1),F(2)],{(0,1):F(-4)})
        self.assertFalse(bad['feasible']);self.assertEqual(F(bad['subset_deficit_Ha']),1)
        with self.assertRaises(ValueError):allocation([F(-1),F(2)],{})
    def test_density_keeps_sixth_degree(self):
        result=closure_demo()['six_mode_allowed_closure']
        self.assertEqual(result['max_degree'],6);self.assertEqual(result['degree6_terms'],5)
        self.assertEqual(F(result['quartic_truncation_L1_operator_penalty']),F(23,147))
    def test_ldl_is_positive_test_not_gap(self):
        self.assertTrue(ldlt_positive([[F(1),F(1,2)],[F(1,2),F(1)]])[0])
        self.assertFalse(ldlt_positive([[F(1),F(1)],[F(1),F(1)]])[0])
        self.assertFalse(ldlt_positive([[F(1),F(2)],[F(2),F(1)]])[0])
    def test_lower_mutations(self):
        v=mono(((1,0),(1,1),(0,3),(0,2)),F(1,7))
        h=canonical(add(*(mono(((1,i),(0,i)),F(e)) for i,e in enumerate([-3,-2,1,2])),v,adj(v)))
        data={'modes':4,'particles':2,'hamiltonian':encode(h)};cert,rec=build(data,2)
        self.assertLessEqual(F(rec['best']['lower_Ha']),F(-5))
        for mutation in ('negative_weight','binding','charge'):
            bad=copy.deepcopy(cert)
            if mutation=='negative_weight':bad['weights'][0]='-1'
            if mutation=='binding':bad['fixture_sha256']='wrong'
            if mutation=='charge':bad['taus'][0]=encode(mono(((1,0),)))
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):check(data,bad)


if __name__=='__main__':unittest.main()
