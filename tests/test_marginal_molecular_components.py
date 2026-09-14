import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from experiments.marginal_molecular_components import components, prepare, schur_matrix, schur_pivots, replay, component_moments
from experiments.marginal_molecular_gershgorin import prepare as scalar_prepare, schur_pivots as scalar_pivots
from experiments.marginal_general_schur import moment_recurrence
from experiments.marginal_symbolic import decode, encode, add, mono

ROOT = Path(__file__).resolve().parents[1]


class MolecularComponentTests(unittest.TestCase):
    def test_stronger_bound_and_uncoupled_low_component_gate(self):
        a = [[F(4), F(0), F(3)], [F(0), F(1), F(0)], [F(3), F(0), F(10)]]
        data = prepare(a, [0])
        self.assertIsNone(scalar_pivots(scalar_prepare(a, [0]), F(0)))
        self.assertEqual(schur_matrix(data, F(0)), [[F(31, 10)]])
        self.assertIsNotNone(schur_pivots(data, F(0)))
        # The uncoupled eigenvalue 1 prevents a global lower bound of 2.
        self.assertIsNone(schur_pivots(data, F(2)))
        self.assertIsNone(schur_pivots(data, F(2), True))

    def test_exact_graph_edges_absolute_radii_and_single_component_equivalence(self):
        tiny = F(1, 10**50)
        self.assertEqual(components([[F(1), tiny], [tiny, F(1)]], [0, 1]), [[0, 1]])
        a = [[F(4), F(1), F(2)], [F(1), F(3), F(1)], [F(2), F(1), F(5)]]
        data = prepare(a, [0])
        self.assertEqual(data['groups'][0]['lower'], F(2))
        self.assertEqual(schur_pivots(data, F(-1)), scalar_pivots(scalar_prepare(a, [0]), F(-1)))
        exact = schur_matrix(data, F(-1), True)[0][0]
        self.assertEqual(exact, F(5) - F(18, 23))
        self.assertGreaterEqual(exact, schur_matrix(data, F(-1))[0][0])

    def test_spectral_thresholds_must_be_verified(self):
        a = [[F(4), F(1), F(2)], [F(1), F(3), F(1)], [F(2), F(1), F(5)]]
        strengthened = prepare(a, [0], ['5/2'])
        self.assertEqual(strengthened['groups'][0]['lower'], F(5, 2))
        for values in ([], ['3'], ['1']):
            with self.assertRaises(ValueError): prepare(a, [0], values)

    def test_reference_partition_covers_Q_and_norm_shift_is_required(self):
        a = [[F(4), F(1), F(2)], [F(1), F(3), F(1)], [F(2), F(1), F(5)]]
        data = prepare(a, [0], partition=[[1], [2]])
        self.assertEqual(data['off_block_norm'], F(1))
        self.assertEqual(schur_matrix(data, F(0), 'moments'), [[F(5, 2)]])
        self.assertIsNone(schur_matrix(data, F(2), 'moments'))
        for partition in ([], [[1]], [[1], [1, 2]], [[0, 1], [2]], [[True], [2]]):
            with self.assertRaises(ValueError): prepare(a, [0], partition=partition)

    def test_connected_stress_replays_with_fixed_small_reference_blocks(self):
        from experiments.marginal_molecular_component_stress import perturb, SOURCE
        source = json.loads(SOURCE.read_text())
        for denominator in (1000, 10000):
            p = ROOT / 'results/marginal_molecular_gershgorin' / ('square_connected_1_' + str(denominator))
            c = json.loads((p / 'certificate.json').read_text())
            self.assertEqual(c['hamiltonian'], perturb(source, F(1, denominator))['hamiltonian'])
            r = replay(c)
            self.assertEqual(F(r['off_block_norm_bound']), F(7, denominator))
            self.assertLessEqual(max(x['degree'] for x in r['response_recurrences']), 7)
            self.assertTrue(all(F(x) > 0 for x in r['positive_denominators']))
            if denominator == 10000: self.assertLess(r['width_float'], 2e-8)
            bad = copy.deepcopy(c)
            bad['reference_components'].pop()
            with self.assertRaises(ValueError): replay(bad)

    def test_moment_response_equals_direct_solve_and_recurrence_is_exact(self):
        a = [[F(4), F(1), F(2)], [F(1), F(3), F(1)], [F(2), F(1), F(5)]]
        data = prepare(a, [0])
        for lower in (F(-2), F(-1), F(0)):
            self.assertEqual(schur_matrix(data, lower, 'moments'), schur_matrix(data, lower, True))
        recurrence = component_moments(data['groups'][0])
        self.assertEqual(recurrence['annihilator'], [F(14), F(-8), F(1)])
        def apply(columns):
            return [{0: 3 * v.get(0, 0) + v.get(1, 0),
                     1: v.get(0, 0) + 5 * v.get(1, 0)} for v in columns]
        with self.assertRaises(ValueError): moment_recurrence([{0: F(1), 1: F(2)}], apply, max_degree=1)
        with self.assertRaises(ValueError): moment_recurrence([], apply, max_degree=True)
        self.assertEqual(moment_recurrence([], apply)['annihilator'], [F(1)])
        uncoupled = prepare([[F(4), F(0)], [F(0), F(1)]], [0])
        self.assertIsNone(schur_matrix(uncoupled, F(2), 'moments'))
        self.assertEqual(schur_matrix(uncoupled, F(0), 'moments'), [[F(4)]])

    def test_moment_replay_does_not_use_inverse_solves(self):
        p = ROOT / 'results/marginal_molecular_gershgorin/square_coupling_components_moments/certificate.json'
        c = json.loads(p.read_text())
        with patch('experiments.marginal_molecular_components.solve_positive', side_effect=AssertionError('Direct inverse path')):
            r = replay(c)
        self.assertLess(r['width_float'], 2e-10)
        self.assertLessEqual(max(x['degree'] for x in r['response_recurrences']), 6)
        bad = copy.deepcopy(c)
        bad['hamiltonian'] = encode(add(decode(c['hamiltonian'], 8, 4), mono((), F(-1))))
        with self.assertRaises(ValueError): replay(bad)

    def test_saved_component_methods_and_false_recipes(self):
        for source in ('square_coupling', 'square_witness'):
            for suffix in ('_components', '_components_spectral', '_components_exact', '_components_moments'):
                p = ROOT / 'results/marginal_molecular_gershgorin' / (source + suffix)
                c = json.loads((p / 'certificate.json').read_text())
                r = replay(c)
                expected = json.loads((p / 'receipt.json').read_text())
                self.assertEqual(r['width'], expected['width'])
                self.assertTrue(all(F(x) > 0 for x in r['positive_denominators']))
                self.assertLessEqual(F(r['width']), F(expected['previous_width']))
        for update in ({'response': 'unknown'}, {'component_bound': 'spectral'},
                       {'component_lowers': []}, {'retained_states': [15, 15]}, {'lower': '10'}):
            bad = copy.deepcopy(c)
            bad.update(update)
            with self.assertRaises(ValueError): replay(bad)


if __name__ == '__main__':
    unittest.main()
