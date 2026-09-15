import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np

from experiments.marginal_implicit_certificate import orthogonal_coordinates, congruence, rayleigh, suggest_upper, suggest_lower, replay, witness_krylov_data, rational_text
from experiments.marginal_schur_transfer import matmul
from experiments.marginal_defect_dicke import enlarged_workspace, workspace_data, targeted_workspace
from experiments.marginal_general_schur import fixture, apply_columns, gram
from experiments.marginal_symbolic import add, decode, encode, mono
from tests.test_marginal_defect_dicke import expand

ROOT = Path(__file__).resolve().parents[1]


class ImplicitCertificateTests(unittest.TestCase):
    def test_large_rational_receipts_preserve_exact_digits_and_budget(self):
        import sys
        old_limit = sys.get_int_max_str_digits()
        self.assertEqual(rational_text(F(10**5000 + 3, 7)), '1' + '0' * 4999 + '3/7')
        self.assertEqual(rational_text(F(-10**5000 - 3, 7)), '-1' + '0' * 4999 + '3/7')
        for value in (F(0), F(-9), F(17, 81), F(1, 10**5000 + 3)):
            if value.denominator.bit_length() < 100:
                self.assertEqual(rational_text(value), str(value))
            else:
                self.assertEqual(rational_text(value), '1/1' + '0' * 4999 + '3')
        with self.assertRaises(ValueError): rational_text(F(1 << 65536))
        self.assertEqual(sys.get_int_max_str_digits(), old_limit)

    def test_exact_transform_and_ill_conditioned_proposal(self):
        l = [[F(1), F(0)], [F(10**20), F(1)]]
        lt = [list(row) for row in zip(*l)]
        g = matmul(l, lt)
        physical = [[F(2), F(1, 3)], [F(1, 3), F(3)]]
        a = matmul(l, matmul(physical, lt))
        t, d = orthogonal_coordinates(g)
        self.assertEqual(d, [1, 1])
        self.assertEqual(congruence(g, t), [[1, 0], [0, 1]])
        self.assertEqual(congruence(a, t), physical)
        data = {'metric': g, 'projected_h': a, 'leakage': [[F(0)] * 2 for _ in range(2)],
                'retained_dimension': 2, 'complement_lower': F(10), 'norm_bound': F(0)}
        coefficients, proposal = suggest_upper(data)
        expected = np.linalg.eigvalsh(np.array(physical, dtype=float))[0]
        self.assertAlmostEqual(float(rayleigh(data, coefficients)), expected, places=12)
        self.assertLessEqual(F(proposal['relative_squared_rounding_error_bound']), F(1, 10**24))
        lower, _ = suggest_lower(data)
        self.assertLess(lower, rayleigh(data, coefficients))
        self.assertLess(float(rayleigh(data, coefficients) - lower), 2e-9)
        huge = dict(data, metric=[[x * 10**800 for x in row] for row in g],
                    projected_h=[[x * 10**800 for x in row] for row in a])
        huge_coefficients, _ = suggest_upper(huge)
        self.assertAlmostEqual(float(rayleigh(huge, huge_coefficients)), expected, places=12)

    def test_bad_metrics_and_variational_coordinates_are_rejected(self):
        for g in ([], [[F(0)]], [[1, 1], [1, 1]], [[1, 1], [0, 1]]):
            with self.assertRaises(ValueError): orthogonal_coordinates(g)
        data = {'metric': [[F(1)]], 'projected_h': [[F(2)]], 'retained_dimension': 1}
        for coordinates in (None, [], [0], [True], ['1'], [1, 2]):
            with self.assertRaises(ValueError): rayleigh(data, coordinates)

    def test_implicit_upper_equals_independent_fock_rayleigh(self):
        h = fixture(F(1, 100))
        workspace = enlarged_workspace(h)
        data = workspace_data(workspace)
        coefficients, _ = suggest_upper(data)
        explicit = {}
        for c, v in zip(coefficients, workspace['basis']):
            for state, amplitude in expand(workspace['model'], v).items():
                explicit[state] = explicit.get(state, F(0)) + c * amplitude
        action = apply_columns(h, [explicit])
        expected = gram([explicit], action)[0][0] / gram([explicit], [explicit])[0][0]
        self.assertEqual(rayleigh(data, coefficients), expected)

    def test_saved_interval_replays_without_explicit_state_helpers(self):
        p = ROOT / 'results/marginal_implicit_certificate/cycle_1_100/certificate.json'
        certificate = json.loads(p.read_text())
        with patch('experiments.marginal_defect_dicke.upper_rayleigh', side_effect=AssertionError('Explicit upper path')):
            with patch('experiments.marginal_general_schur.apply_columns', side_effect=AssertionError('Explicit action path')):
                receipt = replay(certificate)
        self.assertEqual(receipt['stage_dimensions'], [6, 14])
        self.assertLess(receipt['width_float'], 3e-5)
        bad = copy.deepcopy(certificate)
        bad['lower'] = '4'
        with self.assertRaises(ValueError): replay(bad)
        bad = copy.deepcopy(certificate)
        bad['upper_coefficients'] = [0] * 14
        with self.assertRaises(ValueError): replay(bad)
        bad = copy.deepcopy(certificate)
        h = decode(bad['hamiltonian'], 10, 4)
        bad['hamiltonian'] = encode(add(h, mono((), F(-1))))
        with self.assertRaises(ValueError): replay(bad)

    def test_invalid_target_recipe_is_rejected_before_enrichment(self):
        workspace = enlarged_workspace(fixture())
        for coordinates in ([], [0] * 14, [True] * 14):
            with self.assertRaises(ValueError): targeted_workspace(workspace, coordinates)

    def test_saved_targeted_space_is_reconstructed(self):
        path = ROOT / 'results/marginal_implicit_certificate/cycle_1_100_targeted/certificate.json'
        certificate = json.loads(path.read_text())
        receipt = replay(certificate)
        self.assertEqual(receipt['stage_dimensions'], [6, 14, 20])
        self.assertLess(receipt['width_float'], 2e-7)
        certificate['upper_krylov_coefficients'] = [1, 0]
        self.assertEqual(replay(certificate)['upper'], receipt['upper'])

    def test_upper_witness_krylov_moments_match_explicit_actions(self):
        h = fixture(F(1, 100))
        workspace = enlarged_workspace(h)
        coordinates, _ = suggest_upper(workspace_data(workspace))
        small = witness_krylov_data(workspace, coordinates)
        psi = {}
        for coefficient, v in zip(coordinates, workspace['basis']):
            for state, value in expand(workspace['model'], v).items():
                psi[state] = psi.get(state, F(0)) + coefficient * value
        hpsi = apply_columns(h, [psi])[0]
        hhpsi = apply_columns(h, [hpsi])[0]
        self.assertEqual(small['metric'], gram([psi, hpsi], [psi, hpsi]))
        self.assertEqual(small['projected_h'], gram([psi, hpsi], [hpsi, hhpsi]))
        original = rayleigh(workspace_data(workspace), coordinates)
        refined, _ = suggest_upper(small)
        self.assertLessEqual(rayleigh(small, refined), original)
        self.assertNotIn('leakage', small)
        hhhpsi = apply_columns(h, [hhpsi])[0]
        degree_two = witness_krylov_data(workspace, coordinates, steps=2)
        self.assertEqual(degree_two['metric'], gram([psi, hpsi, hhpsi], [psi, hpsi, hhpsi]))
        self.assertEqual(degree_two['projected_h'], gram([psi, hpsi, hhpsi], [hpsi, hhpsi, hhhpsi]))
        with self.assertRaises(ValueError): witness_krylov_data(workspace, coordinates, steps=0)


if __name__ == '__main__':
    unittest.main()
