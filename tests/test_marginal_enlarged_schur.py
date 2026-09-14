import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_enlarged_schur import (closure,solve_positive,projected_action,prepare,pivots,replay,proposed_lower)
from experiments.marginal_general_schur import gram,fixture
from experiments.marginal_symbolic import mono,decode

ROOT=Path(__file__).resolve().parents[1]


class EnlargedSchurTests(unittest.TestCase):
    def test_exact_krylov_closure_detects_rank_and_budget(self):
        h=mono(((1,0),(0,0)));seed=[{1:F(1),2:F(1)}]
        basis=closure(h,seed)
        self.assertEqual(len(basis),2)
        self.assertEqual(len(closure(h,[{1:F(1)},{1:F(2)}])),1)
        with self.assertRaises(ValueError):closure(h,seed,max_dimension=1)
        orthogonal=closure(h,seed,orthogonal=True);g=gram(orthogonal,orthogonal)
        self.assertEqual(g[0][1],0);self.assertEqual(g[1][0],0)
        self.assertTrue(all(max(abs(a) for a in c.values())==1 for c in orthogonal))

    def test_exact_metric_solver_handles_nondiagonal_coordinates(self):
        self.assertEqual(solve_positive([[2,1],[1,2]],[[1],[0]]),[[F(2,3)],[F(-1,3)]])
        self.assertEqual(solve_positive([[2,0],[0,3]],[[1],[1]]),[[F(1,2)],[F(1,3)]])
        with self.assertRaises(ValueError):solve_positive([[1,2],[2,1]],[[1],[1]])
        with self.assertRaises(ValueError):solve_positive([[1,0],[1,1]],[[1],[1]])

    def test_leakage_uses_the_full_metric_projection(self):
        h=mono(((1,0),(0,0)));u=[{1:F(1),2:F(1)}]
        g,a,leak=projected_action(h,u)
        self.assertEqual(g,[[F(2)]]);self.assertEqual(a,[[F(1)]])
        self.assertEqual(leak,[{1:F(1,2),2:F(-1,2)}])
        self.assertEqual(gram(u,leak),[[F(0)]])
        self.assertEqual(gram(leak,leak),[[F(1,2)]])

    def test_scalar_schur_gate_and_numerical_proposal(self):
        data={'metric':[[F(1)]],'projected_h':[[F(1)]],'leakage':[[F(1,100)]],
              'complement_lower':F(2),'norm_bound':F(1,10),'retained_dimension':1}
        self.assertIsNotNone(pivots(data,F(98,100)))
        self.assertIsNone(pivots(data,F(999,1000)))
        proposal=proposed_lower(data)
        self.assertIsNotNone(pivots(data,F(proposal['lower'])))
        self.assertEqual(proposal['acceptance'],'exact rational positivity')

    def test_automatic_enrichment_dimensions_and_remaining_leakage(self):
        data=prepare(fixture(F(1,100)),rounds=2)
        self.assertEqual(data['stage_dimensions'],[6,14,32])
        self.assertEqual(data['retained_support'],232)
        self.assertTrue(any(a for row in data['leakage'] for a in row))
        self.assertTrue(all(not a for row in data['leakage'][:6] for a in row))

    def test_cycle_large_perturbation_certificate(self):
        p=ROOT/'results/marginal_enlarged_schur/cycle_1_100_rounds2/certificate.json'
        certificate=json.loads(p.read_text())
        self.assertEqual(decode(certificate['hamiltonian'],10,4),fixture(F(1,100)))
        r=replay(certificate)
        self.assertGreater(r['width_float'],0);self.assertLess(r['width_float'],1e-7)
        self.assertEqual(r['retained_dimension'],32)

    def test_targeted_mixed_certificate_requires_exact_coefficients(self):
        p=ROOT/'results/marginal_enlarged_schur/mixed_1_100_rounds2_targeted_physical/certificate.json'
        certificate=json.loads(p.read_text())
        self.assertEqual(decode(certificate['hamiltonian'],10,4),fixture(F(1,100),True))
        r=replay(certificate)
        self.assertGreater(r['width_float'],0);self.assertLess(r['width_float'],1e-7)
        self.assertEqual(r['stage_dimensions'],[6,36,42])
        bad=copy.deepcopy(certificate);bad['target_coefficients']=[0]*36
        with self.assertRaises(ValueError):replay(bad)

    def test_false_lower_and_tampered_upper_are_rejected(self):
        p=ROOT/'results/marginal_enlarged_schur/cycle_1_1000/certificate.json'
        certificate=json.loads(p.read_text());bad=copy.deepcopy(certificate);bad['lower']='4'
        with self.assertRaises(ValueError):replay(bad)
        bad=copy.deepcopy(certificate);bad['independent_upper']['amplitudes']=[1]
        with self.assertRaises(ValueError):replay(bad)
        bad=copy.deepcopy(certificate);bad['independent_upper']['numerical_energy']=-1000
        self.assertEqual(replay(bad)['upper'],replay(certificate)['upper'])

    def test_invalid_enrichment_configuration_is_rejected(self):
        h=fixture()
        with self.assertRaises(ValueError):prepare(h,rounds=3)
        with self.assertRaises(ValueError):prepare(h,orthogonal='false')
        with self.assertRaises(ValueError):prepare(h,target_coefficients=[1])


if __name__=='__main__':unittest.main()
