"""Meaningful exact checks at the new spin-changing operator boundary."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import add, canonical, expand_squares, mono, product
from research.joint_patterns_20260913.core import anticommutator
from research.joint_patterns_20260913.spin_diagnostic import spin_generator
from research.molecular_collective_20260913.core import digest, extract, factor_operators
from research.spin_completion_20260913.core import generator, expand_certificate, replay
from research.spin_subspace_20260913.core import combine


class SpinCompletion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]/'results/molecular_collective_20260913/campaign/h4'
        cls.data = json.loads((root/'fixture.json').read_text()); cls.tail = json.loads((root/'rank_6/tail.json').read_text())
        cls.patterns = [q for _, q in factor_operators(extract(cls.data), cls.tail)]
        cls.cert = {'kind': 'spin_completion_subspace_v1', 'fixture_sha256': digest(cls.data), 'tail_sha256': digest(cls.tail),
            'b': '0', 'number_multiplier': [], 'denominator': 7, 'base_blocks': [],
            'anti_blocks': [{'name': 'spin_change', 'generators': [[-1, -1, -1, 0], [0, 0, 1, 1], [1, 1, 0, 0]],
                'directions': [[2, -3, 5], [-1, 4, 2]], 'direction_denominator': 11, 'factor': [[3, -2], [1, 4]]}]}

    def test_diagonal_spin_components_reproduce_previous_operators(self):
        for k in range(len(self.patterns)):
            for spin in (0, 1):
                for mode in range(8):
                    self.assertEqual(generator(self.patterns, [k, spin, spin, mode], 8),
                        spin_generator(self.patterns, [k, spin, mode], 8))

    def test_spin_change_and_its_adjoint_have_correct_CAR_order(self):
        # Explicit off-diagonal spatial density: a0^dag a2 + a2^dag a0,
        # copied in both spins. Its up/down component has a different adjoint.
        q = add(mono(((1, 0), (0, 2))), mono(((1, 2), (0, 0))),
            mono(((1, 1), (0, 3))), mono(((1, 3), (0, 1))))
        B = generator([q], [0, 0, 1, 1], 4)
        expected = product(add(mono(((1, 0), (0, 3))), mono(((1, 2), (0, 1)))), mono(((0, 1),)))
        self.assertEqual(B, expected)
        dagger = product(mono(((1, 1),)), add(mono(((1, 3), (0, 0))), mono(((1, 1), (0, 2)))))
        self.assertEqual(canonical(adj(B)), dagger)
        self.assertLessEqual(max(map(len, anticommutator(B, B))), 4)

    def test_exact_two_stage_factor_matches_direct_mixed_anticommutators(self):
        cert = self.cert; expanded = expand_certificate(self.data, self.tail, cert)
        actual, _ = expand_squares(expanded['blocks'], expanded['denominator'], 8)
        block = cert['anti_blocks'][0]
        polys = [generator(self.patterns, ref, 8) for ref in block['generators']]
        directions = [combine(polys, row, 11) for row in block['directions']]
        expected = {}
        for row in block['factor']:
            q = combine(directions, row, 7); expected = add(expected, anticommutator(q, q))
        self.assertEqual(actual, expected); self.assertLessEqual(max(map(len, actual)), 4)

    def test_tail_transfer_and_input_refusals(self):
        result = replay(self.data, self.tail, self.cert)
        self.assertEqual(F(result['original_lower_Ha']), F(result['retained']['lower'])+F(result['tail']['lower_operator_shift_Ha']))
        for ref in ([0, 2, 0, 1], [-1, 0, 0, 1], [0, 0, 1, 8], [0, 0, True, 1], [0, 0, 1]):
            with self.subTest(ref=ref), self.assertRaises(ValueError):
                generator(self.patterns, ref, 8)
        for key, value in [('directions', [[1, 2]]), ('directions', [[0, 0, 0]]),
                ('direction_denominator', True), ('factor', [[1., 2]])]:
            cert = deepcopy(self.cert); cert['anti_blocks'][0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                expand_certificate(self.data, self.tail, cert)
        cert = deepcopy(self.cert); cert['tail_sha256'] = 'wrong'
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)
        cert = deepcopy(self.cert); cert['base_blocks'] = [{'name': 'invalid', 'words': [[[1, 0], [0, 1], [0, 2]]], 'factor': [[1]]}]
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)


if __name__ == '__main__':
    unittest.main()
