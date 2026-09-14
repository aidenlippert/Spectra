import json
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_adapted_block import effective_model, replay, krylov_lower
from experiments.marginal_boundary_unitary import model_shift, _physical_source
from experiments.marginal_hopping_filter import replay as hopping_replay
from experiments.marginal_symbolic import decode
from experiments.marginal_tiled_upper import _uniform_h
from tests.test_marginal_hopping_filter import action, dot


ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/marginal_graded_hubbard8/tiled_eight_upper/certificate.json'


class AdaptedBlockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads(CERT.read_text())
        cls.c = cls.certificate['upper']
        cls.h = cls.certificate['hamiltonian']

    def test_constant_polynomial_replays_exact_original_upper(self):
        eta = F(57277, 250000)
        adapted = replay(self.c, self.h, eta, [1], 1_000_000)
        hopping = hopping_replay(self.c, self.h, eta, 1_000_000)
        self.assertEqual(adapted['upper'], hopping['upper'])
        self.assertEqual(adapted['degree'], 0)

    def test_nonconstant_recipe_against_direct_determinant_car(self):
        eta=F(1,5)
        oracle,source,_=_physical_source(self.c,self.h)
        full={s:a*phase for r,a in source.items() for s,phase in oracle.orbit(r)[3].items()}
        model,_=effective_model(eta)
        k_full=action(list(decode(model['hamiltonian'],16,4).items()),full)
        # p(x)=1-x/10; x=K+4I. Direct determinant action bypasses
        # both the signed-orbit polynomial recurrence and its energy readout.
        adapted={s:6*full.get(s,0)-k_full.get(s,0) for s in full.keys()|k_full.keys()}
        norm=dot(adapted,adapted)
        h_terms=list(decode(_uniform_h(8)['hamiltonian'],16,4).items())
        e=F(dot(adapted,action(h_terms,adapted)),norm)
        k=F(dot(adapted,action(list(decode(model['hamiltonian'],16,4).items()),adapted)),norm)
        result=replay(self.c,self.h,eta,[1,'-1/10'],24)
        self.assertEqual(F(result['block_energy']),e)
        self.assertEqual(F(result['effective_energy']),k)
        self.assertEqual(F(result['upper']),3*e+2*F(result['merge_shift']))

    def test_effective_symmetry_and_operator_identity_at_zero_and_nonzero_eta(self):
        for eta in (0, F(1, 5)):
            with self.subTest(eta=eta):
                model, weight = effective_model(eta)
                if eta == 0:
                    self.assertEqual(decode(model['hamiltonian'], 16, 4),
                                     decode(_uniform_h(8)['hamiltonian'], 16, 4))
                result = replay(self.c, self.h, eta, [1], 8)
                symmetry = result['effective_symmetry']
                self.assertTrue(symmetry['spin_exchange_exact'])
                self.assertTrue(symmetry['particle_hole_fixed_number_exact'])
                self.assertEqual(result['merge_shift'],
                                 str(model_shift(result['edge_reductions'], eta,
                                             'linear_filter')))
                self.assertEqual(result['thermodynamic_objective_per_block'],
                                 str(F(result['effective_energy']) + 4 * F(weight)
                                     - 2 * F(eta) / (1 + F(eta) * F(eta))))

    def test_degree_zero_krylov_lower_is_exact_source_rayleigh_gate(self):
        accepted = krylov_lower(self.c, self.h, 0, 0, -5)
        self.assertTrue(accepted['accepted'])
        self.assertEqual(accepted['gram_dimension'], 1)
        self.assertTrue(accepted['matrix_psd']['positive_semidefinite'])
        with self.assertRaises(ValueError):
            krylov_lower(self.c, self.h, 0, 0, 0)

    def test_refusal_gates(self):
        with self.assertRaises(ValueError): replay(self.c, self.h, 0.1, [1], 8)
        with self.assertRaises(ValueError): replay(self.c, self.h, 0, [0], 8)
        with self.assertRaises(ValueError): replay(self.c, self.h, 0, [1] + [0] * 9, 8)
        with self.assertRaises(ValueError): replay(self.c, self.h, 0, [1] * 10, 8)
        with self.assertRaises(ValueError): replay(self.c, self.h, 0, [1], 7)
        with self.assertRaises(ValueError): replay(self.c, self.h, 0, [1], 10**9 + 1)


if __name__ == '__main__':
    unittest.main()
