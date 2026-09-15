"""Exact certificate tests at the new two-stage coordinate-map boundary."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

from experiments.marginal_symbolic import add, expand_squares
from research.joint_patterns_20260913.core import anticommutator
from research.joint_patterns_20260913.spin_diagnostic import spin_generator
from research.molecular_collective_20260913.core import digest, extract, factor_operators
from research.spin_subspace_20260913.core import combine, expand_certificate, replay


class SubspaceCertificate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]/'results/molecular_collective_20260913/campaign/h4'
        cls.data = json.loads((root/'fixture.json').read_text()); cls.tail = json.loads((root/'rank_6/tail.json').read_text())
        cls.patterns = [q for _, q in factor_operators(extract(cls.data), cls.tail)]
        cls.cert = {'kind': 'spin_pattern_subspace_v1', 'fixture_sha256': digest(cls.data), 'tail_sha256': digest(cls.tail),
            'b': '0', 'number_multiplier': [], 'denominator': 7, 'base_blocks': [],
            'anti_blocks': [{'name': 'test', 'generators': [[-1, -1, 0], [0, 0, 1], [1, 1, 0]],
                'directions': [[2, -3, 5], [-1, 4, 2]], 'direction_denominator': 11, 'factor': [[3, -2], [1, 4]]}]}

    def test_integer_coordinate_composition_matches_direct_anticommutators(self):
        cert = self.cert; expanded = expand_certificate(self.data, self.tail, cert)
        squares, _ = expand_squares(expanded['blocks'], expanded['denominator'], 8)
        block = cert['anti_blocks'][0]
        polys = [spin_generator(self.patterns, ref, 8) for ref in block['generators']]
        directions = [combine(polys, row, 11) for row in block['directions']]
        expected = {}
        for row in block['factor']:
            q = combine(directions, row, 7); expected = add(expected, anticommutator(q, q))
        self.assertEqual(squares, expected); self.assertLessEqual(max(map(len, squares)), 4)

    def test_same_lower_tail_transfer(self):
        result = replay(self.data, self.tail, self.cert)
        self.assertEqual(F(result['original_lower_Ha']), F(result['retained']['lower'])+F(result['tail']['lower_operator_shift_Ha']))

    def test_bad_directions_factors_and_denominators_refused(self):
        edits = [('directions', [[1, 2]]), ('directions', [[0, 0, 0]]),
            ('directions', [[1, False, 2]]), ('direction_denominator', 0),
            ('direction_denominator', True), ('factor', [[1, 2, 3]]), ('factor', [[1., 2]])]
        for key, value in edits:
            cert = deepcopy(self.cert); cert['anti_blocks'][0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                expand_certificate(self.data, self.tail, cert)

    def test_frozen_binding_and_old_baseline_gates_preserved(self):
        cert = deepcopy(self.cert); cert['tail_sha256'] = 'wrong'
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)
        cert = deepcopy(self.cert); cert['base_blocks'] = [{'name': 'cubic', 'words': [[[1, 0], [0, 1], [0, 2]]], 'factor': [[1]]}]
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)
        cert = deepcopy(self.cert); cert['anti_blocks'][0]['generators'][1] = [0, 2, 1]
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)


if __name__ == '__main__':
    unittest.main()
