"""Independent rational tests for the obstruction's new Gram contraction."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import unittest

from experiments.marginal_symbolic import add, canonical, mono, multiplier_basis, number_shift, product
from experiments.marginal_hunt_car import adj
from research.certificate_scaling.commutator_dual_witness import evaluate, seed
from research.joint_patterns_20260913.core import anticommutator
from research.joint_patterns_20260913.dual import affine_round, integer_grams, check
from research.molecular_collective_20260913.core import digest


def trace_moments(m, n):
    return {tuple((1, i) for i in left)+tuple((0, i) for i in right):
        seed(tuple((1, i) for i in left)+tuple((0, i) for i in right), m, n)
        for k in range(3) for left in combinations(range(m), k) for right in combinations(range(m), k)}


class DualContraction(unittest.TestCase):
    def test_integer_antigram_matches_exact_polynomial_evaluation(self):
        y = trace_moments(4, 2)
        polys = [add(mono(((0, 0),), F(1, 3)), mono(((1, 1), (0, 0), (0, 2)), F(2, 5))),
            add(mono(((0, 3),), F(-2, 7)), mono(((1, 2), (0, 1), (0, 3)), F(3, 11)))]
        matrices, scale = integer_grams(polys, [y], 'anticommutator')
        for i in range(2):
            for j in range(2):
                self.assertEqual(F(matrices[0][i][j], scale), evaluate(anticommutator(polys[i], polys[j]), y))

    def test_integer_square_gram_matches_polynomial_evaluation(self):
        y = trace_moments(4, 2)
        polys = [mono(((1, 0), (0, 1)), F(2, 3)), mono(((1, 2), (0, 2)), F(-1, 5))]
        matrices, scale = integer_grams(polys, [y], 'square')
        for i in range(2):
            for j in range(2):
                self.assertEqual(F(matrices[0][i][j], scale), evaluate(product(canonical(adj(polys[i])), polys[j]), y))

    def test_affine_round_repairs_all_number_identities(self):
        y = trace_moments(4, 2)
        raw = {'rows': list(y), 'values': [float(v)+((i % 7)-3)*1e-7 for i, v in enumerate(y.values())]}
        repaired, _ = affine_round(raw, 4, 2)
        self.assertEqual(repaired[()], 1)
        for q in multiplier_basis(4, max_body=1):
            self.assertEqual(evaluate(product(number_shift(4, 2), q), repaired), 0)
        for w, v in repaired.items():
            self.assertEqual(evaluate(canonical(adj(mono(w))), repaired), v)

    def test_false_normalization_ideal_and_symmetry_claims_refused(self):
        root = Path(__file__).resolve().parents[2]/'results/molecular_collective_20260913/campaign/h6'
        data = json.loads((root/'fixture.json').read_text()); tail = json.loads((root/'rank_10/tail.json').read_text())
        witness = {'kind': 'joint_density_dual_v1', 'fixture_sha256': digest(data),
            'tail_sha256': digest(tail), 'parity_masks': [], 'moments': [{'word': [], 'value': '2'}]}
        with self.assertRaisesRegex(ValueError, 'normalization'):
            check(data, tail, witness)
        witness['moments'][0]['value'] = '1'
        with self.assertRaisesRegex(ValueError, 'number-ideal'):
            check(data, tail, witness)
        witness['parity_masks'] = [1]
        with self.assertRaisesRegex(ValueError, 'parity'):
            check(data, tail, witness)


if __name__ == '__main__':
    unittest.main()
