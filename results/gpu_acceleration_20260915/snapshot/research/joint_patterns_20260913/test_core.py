"""Focused exact regressions for the new pairing and molecular transfer seam."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import add, canonical, expand_squares, mono, product, scale
from experiments.marginal_transfer_verify import apply_word
from research.molecular_collective_20260913.core import digest, extract, factor_operators
from research.joint_patterns_20260913.core import anticommutator, expand_certificate, generator, prepare_anticommutators, replay

ROOT = Path(__file__).resolve().parents[2]


def action(poly, vector):
    out = {}
    for state, value in vector.items():
        for word, coefficient in poly.items():
            image = apply_word(word, state)
            if image is not None:
                target, sign = image
                out[target] = out.get(target, F(0))+coefficient*value*sign
    return {s: v for s, v in out.items() if v}


class JointPatterns(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = ROOT/'results/molecular_collective_20260913/campaign/h4'
        cls.data = json.loads((root/'fixture.json').read_text())
        cls.tail = json.loads((root/'rank_6/tail.json').read_text())
        p = extract(cls.data); cls.patterns = [q for _, q in factor_operators(p, cls.tail)]
        cls.cert = {'kind': 'joint_density_anticommutator_v1',
            'fixture_sha256': digest(cls.data), 'tail_sha256': digest(cls.tail),
            'b': '0', 'number_multiplier': [], 'denominator': 7, 'base_blocks': [],
            'anti_blocks': [{'name': 'test', 'generators': [[-1, 0], [0, 1], [1, 0]],
                'factor': [[2, -3, 1], [0, 1, 4]]}]}

    def test_disjoint_cubic_leading_terms_cancel(self):
        a = mono(((1, 0), (0, 1), (0, 2)))
        b = mono(((1, 3), (0, 4), (0, 5)))
        self.assertTrue(any(len(w) == 6 for w in product(canonical(adj(a)), b)))
        self.assertEqual(anticommutator(a, b), {})

    def test_overlapping_cross_terms_match_independent_fock_actions(self):
        a = add(mono(((1, 0), (0, 1), (0, 2))), mono(((0, 3),), F(2, 3)))
        b = add(mono(((1, 1), (0, 2), (0, 3))), mono(((0, 0),), F(-1, 5)))
        q = add(scale(a, F(3, 7)), scale(b, F(-2, 11)))
        vector = {s: F((s % 5)-2, 13) for s in range(16)}
        anti = anticommutator(q, q)
        lhs = sum(vector[s]*v for s, v in action(anti, vector).items())
        rhs = sum(v*v for v in action(q, vector).values())+sum(v*v for v in action(canonical(adj(q)), vector).values())
        self.assertEqual(lhs, rhs); self.assertGreaterEqual(rhs, 0)

    def test_even_generators_refused(self):
        with self.assertRaises(ValueError):
            anticommutator(mono(((1, 0), (0, 0))), mono(((0, 1),)))
        with self.assertRaises(ValueError):
            prepare_anticommutators([{'polynomials': [mono(())]}])

    def test_compact_pair_expands_to_exact_adjoint_squares(self):
        expanded = expand_certificate(self.data, self.tail, self.cert)
        squares, _ = expand_squares(expanded['blocks'], expanded['denominator'], 8)
        expected = {}
        block = self.cert['anti_blocks'][0]
        polys = [generator(self.patterns, ref, 8) for ref in block['generators']]
        for row in block['factor']:
            q = add(*(scale(p, F(c, 7)) for p, c in zip(polys, row)))
            expected = add(expected, anticommutator(q, q))
        self.assertEqual(squares, expected)
        self.assertLessEqual(max(map(len, squares)), 4)

    def test_transfer_uses_lower_tail_endpoint(self):
        result = replay(self.data, self.tail, self.cert)
        self.assertEqual(F(result['original_lower_Ha']), F(result['retained']['lower'])+F(result['tail']['lower_operator_shift_Ha']))
        self.assertEqual(result['many_body_states_enumerated'], 0)

    def test_mismatched_tail_and_fixture_refused(self):
        for key in ('tail_sha256', 'fixture_sha256'):
            cert = deepcopy(self.cert); cert[key] = 'bad'
            with self.assertRaises(ValueError):
                expand_certificate(self.data, self.tail, cert)

    def test_invalid_generator_and_factor_refused(self):
        for ref in ([6, 0], [0, 8], [True, 0], [-2, 0]):
            cert = deepcopy(self.cert); cert['anti_blocks'][0]['generators'][0] = ref
            with self.assertRaises(ValueError):
                expand_certificate(self.data, self.tail, cert)
        cert = deepcopy(self.cert); cert['anti_blocks'][0]['factor'][0][0] = 0.5
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)

    def test_baseline_degree_budget_enforced(self):
        cert = deepcopy(self.cert)
        cert['base_blocks'] = [{'name': 'too_large', 'words': [((1, 0), (0, 1), (0, 2))], 'factor': [[1]]}]
        with self.assertRaises(ValueError):
            expand_certificate(self.data, self.tail, cert)


if __name__ == '__main__':
    unittest.main()
