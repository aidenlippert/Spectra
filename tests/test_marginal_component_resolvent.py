import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from experiments.marginal_component_resolvent import polynomial_inverse, full_recurrence
from experiments.marginal_enlarged_schur import solve_positive
from experiments.marginal_molecular_components import prepare, schur_matrix, replay, strongest_edge_partition
from experiments.marginal_schur_transfer import ldl_pivots

ROOT = Path(__file__).resolve().parents[1]


class ComponentResolventTests(unittest.TestCase):
    def test_noncommuting_response_is_upper_bound(self):
        # P has dimension two; Q has a 2x2 block and a singleton.
        a = [[F(x) for x in row] for row in (
            (8, 1, 1, 2, 3), (1, 9, 2, -1, 1),
            (1, 2, 6, 1, 1), (2, -1, 1, 7, -1), (3, 1, 1, -1, 9))]
        data = prepare(a, [0, 1], partition=[[2, 3], [4]])
        actual = prepare(a, [0, 1])
        for b in (F(-2), F(0), F(1)):
            exact = schur_matrix(actual, b, True)
            bounded = schur_matrix(data, b, 'second_order')
            difference = [[x - y for x, y in zip(row, bounded[i])] for i, row in enumerate(exact)]
            self.assertIsNotNone(ldl_pivots(difference))

    def test_zero_R_recovers_exact_response_and_uncoupled_gate(self):
        a = [[F(x) for x in row] for row in ((5, 1, 2, 0), (1, 6, 1, 0), (2, 1, 7, 0), (0, 0, 0, 1))]
        data = prepare(a, [0])
        self.assertEqual(schur_matrix(data, F(0), 'second_order'), schur_matrix(data, F(0), True))
        for b in (F(1), F(2)):
            self.assertIsNone(schur_matrix(data, b, 'second_order'))
        connected = prepare([[F(x) for x in row] for row in ((5, 1, 0), (1, 4, 1), (0, 1, 3))], [0], partition=[[1], [2]])
        self.assertIsNone(schur_matrix(connected, F(2), 'second_order'))

    def test_full_polynomial_inverse_and_false_polynomial(self):
        group = {'matrix': [[F(3), F(1)], [F(1), F(5)]]}
        self.assertEqual(full_recurrence(group)['annihilator'], [F(14), F(-8), F(1)])
        for z in (F(-1), F(0), F(1)):
            shifted = [[x - z * (i == j) for j, x in enumerate(row)] for i, row in enumerate(group['matrix'])]
            self.assertEqual(polynomial_inverse(group, z), solve_positive(shifted, [[F(1), F(0)], [F(0), F(1)]]))
        group['full_recurrence']['annihilator'][0] += 1
        with self.assertRaises(ValueError): polynomial_inverse(group, F(0))
        with self.assertRaises(ValueError): polynomial_inverse({'matrix': [[F(2)]]}, F(2))

    def test_saved_second_order_replay_and_tampering(self):
        for denominator in (1000, 10000):
            directory = ROOT / 'results/marginal_molecular_gershgorin' / f'square_connected_1_{denominator}_second_order'
            c = json.loads((directory / 'certificate.json').read_text())
            with patch('experiments.marginal_molecular_components.solve_positive', side_effect=AssertionError('Direct inverse')):
                r = replay(c)
            saved = json.loads((directory / 'receipt.json').read_text())
            self.assertEqual(r['width'], saved['width'])
            self.assertLess(r['width_float'], 1e-7)
            self.assertLess(r['width_float'], saved['previous_width_float'])
            self.assertLessEqual(max(x['degree'] for x in r['response_recurrences']), 7)
            for update in ({'lower': '10'}, {'reference_components': c['reference_components'][:-1]}):
                bad = copy.deepcopy(c)
                bad.update(update)
                with self.assertRaises(ValueError): replay(bad)

    def test_automatic_partition_uses_actual_edges_and_enforces_cap(self):
        a = [[F(x) for x in row] for row in ((5, 9, 1, 0), (9, 5, 2, 1), (1, 2, 5, 8), (0, 1, 8, 5))]
        self.assertEqual(strongest_edge_partition(a, list(range(4)), 2), [[0, 1], [2, 3]])
        self.assertEqual(strongest_edge_partition(a, list(reversed(range(4))), 2), [[0, 1], [2, 3]])
        for cap in (0, True, 13, 1.5):
            with self.assertRaises(ValueError): strongest_edge_partition(a, list(range(4)), cap)
        for denominator in (1000, 10000):
            base = ROOT / 'results/marginal_molecular_gershgorin' / f'square_connected_1_{denominator}'
            directory = base.with_name(base.name + '_automatic_second_order')
            c = json.loads((directory / 'certificate.json').read_text())
            original = json.loads((base / 'certificate.json').read_text())
            self.assertEqual(c['hamiltonian'], original['hamiltonian'])
            self.assertEqual(c['independent_upper'], original['independent_upper'])
            with patch('experiments.marginal_molecular_components.solve_positive', side_effect=AssertionError('Direct inverse')):
                r = replay(c)
            self.assertLess(r['width_float'], 3e-10)
            self.assertLessEqual(max(r['component_dimensions']), 7)
            self.assertEqual(r['width'], json.loads((directory / 'receipt.json').read_text())['width'])
            for recipe in ({'method': 'strongest_edges', 'max_size': 1}, {'method': 'other', 'max_size': 7}):
                bad = copy.deepcopy(c)
                bad['reference_selection'] = recipe
                with self.assertRaises(ValueError): replay(bad)


if __name__ == '__main__':
    unittest.main()
