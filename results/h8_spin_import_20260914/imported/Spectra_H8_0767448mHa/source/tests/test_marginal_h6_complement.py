import copy
import contextlib
from fractions import Fraction as F
import json
import io
from pathlib import Path
import unittest
import tempfile
from unittest.mock import patch

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_fixed_point_ldl import propose, verify, rounded_division
from experiments.marginal_h6_complement import replay_comparison, replay_factor, checked_blocks
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_sparse_response import replay as response_replay, replay_rank_obstruction, refine_factor, refine_response
from experiments.marginal_symbolic import add, mono, encode

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/marginal_h6'


def fixture(a):
    h = add(*(mono(((1, i), (0, j)), F(x)) for i, row in enumerate(a) for j, x in enumerate(row) if x))
    return {'modes': len(a), 'particles': 1, 'hamiltonian': encode(h), 'retained_states': [1]}


def toy():
    a = [[0, F(1, 10), F(1, 10), F(1, 10)],
         [F(1, 10), 2, 1, 1], [F(1, 10), 1, 2, 1], [F(1, 10), 1, 1, 2]]
    c = fixture(a)
    gamma = F(1, 2)
    q = [[F(x) - gamma * (i == j) for j, x in enumerate(row[1:])] for i, row in enumerate(a[1:])]
    c.update(kind='factor_complement_bound_v1', complement_lower=str(gamma),
             blocks=[{'states': [2, 4, 8], 'factor': propose(q, 10**12)}])
    return a, c


class FixedPointFactorTests(unittest.TestCase):
    def test_signed_nearest_even_rounding(self):
        self.assertEqual([rounded_division(i, 2) for i in range(-5, 6)], [-2, -2, -2, -1, 0, 0, 0, 1, 2, 2, 2])
        with self.assertRaises(ValueError): rounded_division(1, 0)

    def test_independent_residual_and_inverse_norm_arithmetic(self):
        a = [[F(x) for x in row] for row in ((5, -2, 1), (-2, 4, -1), (1, -1, 3))]
        factor = propose(a, 1000)
        r = verify(a, factor)
        n, s = len(a), factor['scale']
        l = [[F(factor['lower'][i][j], s) if j < i else F(i == j) for j in range(n)] for i in range(n)]
        d = [F(x, s) for x in factor['diagonal']]
        residual = [[a[i][j] - sum(l[i][k] * d[k] * l[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
        self.assertEqual(r['residual_norm_bound'], max(sum(abs(x) for x in row) for row in residual))
        inverse = [[F(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                inverse[i][j] = F(i == j) - sum(l[i][k] * inverse[k][j] for k in range(i))
        actual_norm_product = max(sum(abs(x) for x in row) for row in inverse) * max(sum(abs(inverse[i][j]) for i in range(n)) for j in range(n))
        self.assertGreaterEqual(r['inverse_norm_squared_bound'], actual_norm_product)
        self.assertIsNotNone(ldl_pivots([[x - r['margin'] / 2 * (i == j) for j, x in enumerate(row)] for i, row in enumerate(a)]))

    def test_false_factors_nonpositive_matrices_and_invalid_scales(self):
        a = [[F(2), F(1)], [F(1), F(2)]]
        good = propose(a, 1000)
        for update in ({'scale': True}, {'scale': 0}, {'scale': 10**65}, {'diagonal': [2000, 0]},
                       {'diagonal': [200000, 200000]}, {'lower': [[], [1000000]]}, {'lower': [[], [1.0]]}):
            with self.assertRaises(ValueError): verify(a, dict(good, **update))
        for bad in ([[1, 2], [2, 1]], [[0, 0], [0, 1]]):
            with self.assertRaises(ValueError): propose(bad)
        with self.assertRaises(ValueError): verify([[F(1, 3)]], {'scale': 10, 'lower': [[]], 'diagonal': [3]})
        with self.assertRaises(ValueError): propose([[1, 1], [0, 1]])


class ComplementProofTests(unittest.TestCase):
    def test_opt_in_large_block_preserves_coverage_and_physical_gates(self):
        from experiments.marginal_spin_reduction import SpinZeroOracle
        from experiments.marginal_spin_constructor import spin_states
        oracle=SpinZeroOracle({'modes':12,'particles':6,'hamiltonian':[]})
        states=spin_states(oracle);p=states[:32];q=states[32:]
        blocks=[{'states':q}]
        with self.assertRaises(ValueError): checked_blocks(oracle,p,blocks)
        a=checked_blocks(oracle,p,blocks,max_block_dimension=384)
        self.assertEqual(len(a[0]),368)
        for bad in ([{'states':q[:-1]}],[{'states':[p[0]]+q[1:]}],[{'states':q[:-1]+[q[0]]}]):
            with self.assertRaises(ValueError): checked_blocks(oracle,p,bad,max_block_dimension=384)
        for bound in (True,385):
            with self.assertRaises(ValueError): checked_blocks(oracle,p,blocks,max_block_dimension=bound)

    def test_frustrated_triangle_comparison_obstruction_with_physical_gap(self):
        _, c = toy()
        self.assertEqual(replay_factor(c)['complement_lower'], '1/2')
        c.update(kind='comparison_gershgorin_obstruction_v1', target_lower='1/2',
                 comparison_witness={'states': [2, 4, 8], 'amplitudes': [1, 1, 1]})
        r = replay_comparison(c)
        self.assertEqual(F(r['comparison_rayleigh']), 0)
        for update in ({'target_lower': '0'}, {'retained_states': [1, 1]},
                       {'comparison_witness': {'states': [1], 'amplitudes': [1]}},
                       {'comparison_witness': {'states': [2], 'amplitudes': [0]}}):
            with self.assertRaises(ValueError): replay_comparison(dict(c, **update))

    def test_complete_coverage_and_interblock_couplings(self):
        _, c = toy()
        oracle = DeterminantOracle(c)
        for groups in ([[2, 4]], [[2, 4, 4]], [[1, 2, 4]], [[2], [4, 8]]):
            with self.assertRaises(ValueError): checked_blocks(oracle, [1], [{'states': x} for x in groups])
        a = [[0, 0, 0], [0, 2, 0], [0, 0, -1]]
        bad = fixture(a)
        bad.update(kind='factor_complement_bound_v1', complement_lower='0',
                   blocks=[{'states': [2], 'factor': propose([[2]])}])
        with self.assertRaises(ValueError): replay_factor(bad)
        bad['blocks'].append({'states': [4], 'factor': propose([[1]])})
        with self.assertRaises(ValueError): replay_factor(bad)

    def test_factor_gap_and_response_use_the_same_physical_operator(self):
        a, c = toy()
        c.update(kind='factor_response_schur_v1', lower='-1/100',
                 response_basis=[{'states': [s], 'amplitudes': [1]} for s in (2, 4, 8)],
                 independent_upper={'states': [1], 'amplitudes': [1]})
        self.assertEqual(F(response_replay(c)['width']), F(1, 100))
        self.assertIsNotNone(ldl_pivots([[F(x) + F(1, 100) * (i == j) for j, x in enumerate(row)] for i, row in enumerate(a)]))
        for update in ({'lower': '0'}, {'complement_lower': '2'}, {'blocks': []}, {'response_basis': []}):
            with self.assertRaises(ValueError): response_replay(dict(c, **update))
        a[0][0] = -5
        with self.assertRaises(ValueError): response_replay(dict(c, hamiltonian=fixture(a)['hamiltonian']))

    def test_saved_h6_proofs_are_bound_and_replay_without_proposals(self):
        path = DATA / 'complement_proofs'
        fixture_data = json.loads((DATA / 'fixture.json').read_text())
        for name, replay in (('comparison', replay_comparison), ('factor', replay_factor)):
            c = json.loads((path / (name + '_certificate.json')).read_text())
            for key in ('modes', 'particles', 'hamiltonian'):
                self.assertEqual(c[key], fixture_data[key])
            with patch('experiments.marginal_h6_complement.propose', side_effect=AssertionError('Proposal during replay')):
                r = replay(c)
            self.assertEqual(r, json.loads((path / (name + '_receipt.json')).read_text()))
        self.assertEqual(sum(r['block_dimensions']), 892)
        self.assertEqual(max(r['block_dimensions']), 200)
        bad = copy.deepcopy(c)
        bad['blocks'][0]['factor']['diagonal'][0] = 0
        with self.assertRaises(ValueError): replay_factor(bad)

    def test_response_rank_obstruction_and_false_targets(self):
        _, c = toy()
        c.update(kind='response_rank_obstruction_v1', target_lower='-1/100', negative_principal_indices=[0])
        self.assertEqual(replay_rank_obstruction(c)['minimum_response_dimension'], 1)
        for update in ({'target_lower': '-10'}, {'target_lower': '1'},
                       {'negative_principal_indices': [0, 0]}, {'negative_principal_indices': [True]},
                       {'negative_principal_indices': []}, {'negative_principal_indices': [1]}):
            with self.assertRaises(ValueError): replay_rank_obstruction(dict(c, **update))
        c = json.loads((DATA / 'response_rank_obstruction/certificate.json').read_text())
        r = replay_rank_obstruction(c)
        self.assertEqual(r['minimum_response_dimension'], 31)
        self.assertEqual(r, json.loads((DATA / 'response_rank_obstruction/receipt.json').read_text()))

    def test_factor_response_proposal_produces_checked_small_interval(self):
        _, c = toy()
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'factor'
            source.mkdir()
            path, upper = source / 'certificate.json', source / 'upper.json'
            path.write_text(json.dumps(c))
            upper.write_text(json.dumps({'independent_upper': {'states': [1, 2, 4, 8], 'amplitudes': [400, -10, -10, -10]}}))
            with contextlib.redirect_stdout(io.StringIO()): refine_factor(path, upper)
            result = source.with_name('factor_response')
            got = json.loads((result / 'certificate.json').read_text())
            r = response_replay(got)
            self.assertLess(r['width_float'], 1e-7)
            self.assertEqual(r['response_dimension'], 1)
            self.assertEqual(r['unique_action_states'], 4)
            resumed = Path(directory) / 'resumed'
            with contextlib.redirect_stdout(io.StringIO()):
                refine_response(got, resumed, r, F(1, 10**7), str(result))
            self.assertEqual(json.loads((resumed / 'receipt.json').read_text())['width'], r['width'])

    def test_saved_h6_global_interval_and_full_coverage_accounting(self):
        path = DATA / 'complement_proofs_response_continued'
        c = json.loads((path / 'certificate.json').read_text())
        def bounded_solve(a, b):
            from experiments.marginal_enlarged_schur import solve_positive
            self.assertLessEqual(len(a), 32)
            return solve_positive(a, b)
        with patch('experiments.marginal_sparse_response.solve_positive', side_effect=bounded_solve), \
             patch('experiments.marginal_sparse_response.approximate_response', side_effect=AssertionError('Proposal during replay')):
            r = response_replay(c)
        self.assertEqual(r['width'], json.loads((path / 'receipt.json').read_text())['width'])
        self.assertLess(r['width_float'], 2e-10)
        self.assertEqual(r['response_dimension'], 32)
        self.assertEqual(r['unique_action_states'], 924)
        self.assertEqual(r['referenced_determinants'], 924)
        reference = F(json.loads((DATA / 'reference_upper.json').read_text())['upper'])
        self.assertLess(F(r['lower']), reference)
        self.assertLessEqual(reference, F(r['upper']))
        with self.assertRaises(ValueError): response_replay(dict(c, kind='sparse_response_schur_v1'))


if __name__ == '__main__':
    unittest.main()
