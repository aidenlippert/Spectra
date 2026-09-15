from fractions import Fraction as F
import copy
import json
from pathlib import Path
import unittest

from experiments.marginal_general_schur import fixture
from experiments.marginal_targeting_bound import parameters, reference_excitation, internal_residual, check_step, propose_retained_lower


class TargetingBoundTests(unittest.TestCase):
    def test_exact_step_gate_accepts_and_rejects_inaccurate_vectors(self):
        data = {'metric': [[F(1), F(0)], [F(0), F(1)]],
                'projected_h': [[F(0), F(0)], [F(0), F(2)]], 'retained_dimension': 2}
        guarantee = {'witness_upper': '1', 'spectral_width_bound': '4',
                     'contraction_factor': '7/8', 'error_floor': '1/100'}
        accepted = check_step(data, [1, 0], data, [1, 0], F(-1, 10000), guarantee)
        self.assertLessEqual(F(accepted['additive_error_upper']), F(accepted['allowed_additive_error']))
        self.assertEqual(internal_residual(data, [1, 1]), (F(1), F(1)))
        self.assertEqual(internal_residual(data, [17, 17]), (F(1), F(1)))
        for current, following, lower in (([1, 1], [1, 0], F(-1, 10000)),
                                           ([1, 0], [0, 1], F(-1, 10000)),
                                           ([1, 0], [1, 0], F(1, 10000)),
                                           ([1, 0], [1, 0], F(-1))):
            with self.assertRaises(ValueError): check_step(data, current, data, following, lower, guarantee)
        lower = propose_retained_lower(data, [1, 0], guarantee)
        check_step(data, [1, 0], data, [1, 0], lower, guarantee)
        tiny = dict(guarantee, error_floor='1/1000000000000000000000000000000')
        with self.assertRaises(ValueError): propose_retained_lower(data, [1, 0], tiny)

    def test_saved_chain_convergence_is_recomputed_and_tampering_rejected(self):
        from experiments.marginal_implicit_certificate import certificate_workspace, replay
        root = Path(__file__).resolve().parents[1]
        c = json.loads((root / 'results/marginal_implicit_certificate/cycle_1_100_targeted_compressed/certificate.json').read_text())
        w, _, _ = certificate_workspace(c, certify_floor=F(1, 10**7))
        c['targeting_guarantee'] = w['targeting_guarantee']
        result = replay(c)
        self.assertEqual(result['targeting_convergence'], w['targeting_receipt'])
        self.assertEqual(len(result['targeting_convergence']['steps']), 1)
        for modification in ({'retained_lowers': []}, {'upper': '0'}, {'retained_lowers': ['100']},
                             {'error_floor': '0'}, {'error_floor': '1/1000000000000000000000000000000'}):
            bad = copy.deepcopy(c)
            bad['targeting_guarantee'].update(modification)
            with self.assertRaises(ValueError): replay(bad)
        bad = copy.deepcopy(c)
        bad['upper_krylov_coefficients'] = [1, 1]
        with self.assertRaises(ValueError): replay(bad)
        c.pop('targeting_guarantee')
        with self.assertRaises(ValueError): certificate_workspace(c, certify_floor=F(0))

    def test_middle_target_is_bound_to_its_actual_reconstructed_vector(self):
        from experiments.marginal_implicit_certificate import replay
        root = Path(__file__).resolve().parents[1]
        p = root / 'results/marginal_implicit_certificate/cycle_1_100_targeted_convergence_gated_enriched/certificate.json'
        c = json.loads(p.read_text())
        self.assertEqual(len(c['target_chain']), 2)
        c['target_chain'][1] = [1] + [0] * (len(c['target_chain'][1]) - 1)
        with self.assertRaises(ValueError): replay(c)

    def test_excitation_inertia_and_open_gap_rejection(self):
        for pairs in range(2, 9):
            c, pivots = reference_excitation(pairs)
            self.assertEqual(sum(x < 0 for x in pivots), 1)
            self.assertEqual(len(pivots), pairs + 1)
        for pairs, upper in ((5, F(3281, 1000)), (6, F(5375, 1000))):
            for mixed in (False, True):
                h = fixture(F(1, 100), mixed, 2 * pairs)
                p = parameters(h, pairs, upper)
                self.assertTrue(0 < F(p['contraction_factor']) < 1)
                with self.assertRaises(ValueError):
                    parameters(h, pairs, F(p['excited_lower']))

    def test_residual_step_and_variance_with_internal_error(self):
        # Rational diagonal Hamiltonian; retained space need not be spectral.
        # v=(3,4,0)/5 and w=(4,-3,5)/sqrt(50) span U.
        h = [F(0), F(2), F(7)]
        v = [F(3, 5), F(4, 5), F(0)]
        w = [F(4), F(-3), F(5)]
        dot = lambda a, b: sum(x * y for x, y in zip(a, b))
        hv = [x * y for x, y in zip(h, v)]
        mu = dot(v, hv)
        t = [x * dot(w, hv) / dot(w, w) for x in w]
        s = [a - mu * b - d for a, b, d in zip(hv, v, t)]
        self.assertEqual(dot(v, s), 0)
        self.assertEqual(dot(w, s), 0)
        variance = dot(s, s) + dot(t, t)
        self.assertGreaterEqual(variance, mu * (2 - mu))
        trial = [a - b / 7 for a, b in zip(v, s)]
        energy = sum(x * y * y for x, y in zip(h, trial)) / dot(trial, trial)
        self.assertLessEqual(energy, mu - dot(s, s) / 14)
        factor = 1 - (2 - mu) / 14
        self.assertLessEqual(energy, factor * mu + dot(t, t) / 14)


if __name__ == '__main__':
    unittest.main()
