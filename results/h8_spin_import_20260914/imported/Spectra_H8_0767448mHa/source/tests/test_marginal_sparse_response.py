import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_component_resolvent import multiply, transpose
from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_enlarged_schur import solve_positive
from experiments.marginal_sparse_response import prepare, response_matrix, parse_basis, replay, refine
from experiments.marginal_symbolic import add, encode, mono

ROOT = Path(__file__).resolve().parents[1]


def oracle_for(a):
    h = add(*(mono(((1, i), (0, j)), F(x)) for i, row in enumerate(a) for j, x in enumerate(row) if x))
    return DeterminantOracle({'modes': len(a), 'particles': 1, 'hamiltonian': encode(h)})


class SparseResponseTests(unittest.TestCase):
    def assert_psd2(self, a):
        self.assertEqual(a[0][1], a[1][0])
        self.assertGreaterEqual(a[0][0], 0)
        self.assertGreaterEqual(a[1][1], 0)
        self.assertGreaterEqual(a[0][0] * a[1][1] - a[0][1] ** 2, 0)

    def test_noncommuting_bound_monotonicity_and_exact_full_response(self):
        a = [[F(x) for x in row] for row in (
            (8, 1, 1, 2, 3), (1, 9, 2, -1, 1),
            (1, 2, 6, 1, 1), (2, -1, 1, 7, -1), (3, 1, 1, -1, 9))]
        oracle = oracle_for(a)
        q = [row[2:] for row in a[2:]]
        w = [row[:2] for row in a[2:]]
        for b in (F(-1), F(0), F(1)):
            shifted = [[x - b * (i == j) for j, x in enumerate(row)] for i, row in enumerate(q)]
            response = multiply(transpose(w), solve_positive(shifted, w))
            actual = [[a[i][j] - b * (i == j) - response[i][j] for j in range(2)] for i in range(2)]
            previous = response_matrix(prepare(oracle, [1, 2], F(3), []), b)
            basis = []
            for state in (4, 8, 16):
                basis.append({state: F(1)})
                current = response_matrix(prepare(oracle, [1, 2], F(3), basis), b)
                self.assert_psd2([[current[i][j] - previous[i][j] for j in range(2)] for i in range(2)])
                self.assert_psd2([[actual[i][j] - current[i][j] for j in range(2)] for i in range(2)])
                previous = current
            self.assertEqual(previous, actual)

    def test_optimized_expression_equals_explicit_residual_bound(self):
        a = [[F(x) for x in row] for row in ((8, 1, 2), (1, 5, 1), (2, 1, 7))]
        oracle = oracle_for(a)
        data = prepare(oracle, [1], F(3), [{2: F(1), 4: F(2)}])
        b, delta = F(1), F(2)
        shifted = [[F(4), F(1)], [F(1), F(6)]]
        k, w = [[F(1)], [F(2)]], [[F(1)], [F(2)]]
        ak = multiply(shifted, k)
        d = multiply(transpose(ak), ak)[0][0] - delta * multiply(transpose(k), ak)[0][0]
        l = multiply(transpose(ak), w)[0][0] - delta * multiply(transpose(k), w)[0][0]
        x = [[row[0] * l / d] for row in k]
        ax = multiply(shifted, x)
        residual = [[w[i][0] - ax[i][0]] for i in range(2)]
        explicit = (2 * multiply(transpose(w), x)[0][0] - multiply(transpose(x), ax)[0][0]
                    + multiply(transpose(residual), residual)[0][0] / delta)
        self.assertEqual(response_matrix(data, b), [[F(8) - b - explicit]])

    def test_gap_rank_and_response_support_gates(self):
        oracle = oracle_for([[5, 1], [1, 2]])
        data = prepare(oracle, [1], F(2), [{2: F(1)}])
        self.assertIsNone(response_matrix(data, F(2)))
        # Independent K, but A(A-delta I) has zero curvature: no inverse allowed.
        with self.assertRaises(ValueError): response_matrix(data, F(0))
        with self.assertRaises(ValueError): prepare(oracle, [1], F(1), [{2: F(1)}, {2: F(2)}])
        for recipe in ([], [{'states': [1], 'amplitudes': [1]}],
                       [{'states': [2, 2], 'amplitudes': [1, 1]}],
                       [{'states': [2], 'amplitudes': [0]}],
                       [{'states': [2], 'amplitudes': [1.]}]):
            with self.assertRaises(ValueError): parse_basis(oracle, {1}, recipe)

    def test_saved_replay_uses_only_small_response_solve(self):
        p = ROOT / 'results/marginal_sparse_molecular/square_residual_response'
        c = json.loads((p / 'certificate.json').read_text())
        def small_solve(a, b):
            self.assertLessEqual(len(a), 16)
            return solve_positive(a, b)
        with patch('experiments.marginal_sparse_response.solve_positive', side_effect=small_solve), \
             patch('experiments.marginal_sparse_response.approximate_response', side_effect=AssertionError('Proposal during replay')), \
             patch('experiments.marginal_molecular_gershgorin.matrix', side_effect=AssertionError('Full matrix')):
            r = replay(c)
        self.assertLess(r['width_float'], 2e-10)
        self.assertEqual(r['response_dimension'], 12)
        self.assertEqual(r['unique_action_states'], 70)
        self.assertEqual(r['width'], json.loads((p / 'receipt.json').read_text())['width'])
        for update in ({'response_basis': []}, {'lower': '10'}, {'complement_tree': ['bound']},
                       {'response_basis': [c['response_basis'][0], c['response_basis'][0]]}):
            bad = copy.deepcopy(c)
            bad.update(update)
            with self.assertRaises(ValueError): replay(bad)

    def test_proposer_reconstructs_without_full_sector_eigensolve(self):
        import numpy as np
        c = json.loads((ROOT / 'results/marginal_sparse_molecular/square_residual/certificate.json').read_text())
        original_eigh, original_solve = np.linalg.eigh, np.linalg.solve
        def small_eigh(a):
            self.assertLessEqual(len(a), 32)
            return original_eigh(a)
        def small_solve(a, b):
            self.assertLessEqual(len(a), 16)
            return original_solve(a, b)
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source'
            source.mkdir()
            (source / 'certificate.json').write_text(json.dumps(c))
            with patch('numpy.linalg.eigh', side_effect=small_eigh), patch('numpy.linalg.solve', side_effect=small_solve), \
                 patch('experiments.marginal_molecular_gershgorin.matrix', side_effect=AssertionError('Full matrix')), \
                 contextlib.redirect_stdout(io.StringIO()):
                refine(source / 'certificate.json')
            receipt = json.loads((source.with_name('source_response') / 'receipt.json').read_text())
            self.assertLess(receipt['width_float'], 2e-10)
            self.assertEqual(receipt['construction_unique_action_states'], 70)
            self.assertGreaterEqual(receipt['construction_unique_action_states'], receipt['proposal_unique_action_states'])


if __name__ == '__main__':
    unittest.main()
