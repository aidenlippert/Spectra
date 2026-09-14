import copy
import contextlib
from fractions import Fraction as F
import json
import io
from math import comb
from pathlib import Path
import unittest
import tempfile
from unittest.mock import patch

from experiments.marginal_determinant_tree import DeterminantOracle, GapFailure
from experiments.marginal_sparse_molecular import replay
from experiments.marginal_molecular_gershgorin import matrix, prepare, SOURCES
from experiments.marginal_symbolic import encode, mono, add

ROOT = Path(__file__).resolve().parents[1]


class DeterminantTreeTests(unittest.TestCase):
    def test_actions_diagonal_and_every_prefix_bound_against_full_oracle(self):
        for name, path in SOURCES.items():
            c = json.loads((ROOT / path).read_text())
            states, a = matrix(c)
            oracle = DeterminantOracle(c)
            retained = set(states[::3])
            q = [i for i, s in enumerate(states) if s not in retained]
            rows = {states[i]: a[i][i] - sum(abs(a[i][j]) for j in q if j != i) for i in q}
            for j, state in enumerate(states):
                self.assertEqual(oracle.action(state), {s: a[i][j] for i, s in enumerate(states) if a[i][j]})
                self.assertEqual(oracle.diagonal_lower(255, state), a[j][j])
            for depth in range(9):
                mask = (1 << depth) - 1
                for bits in range(1 << depth):
                    completions = [i for i, s in enumerate(states) if s & mask == bits]
                    if not completions:
                        continue
                    self.assertLessEqual(oracle.diagonal_lower(mask, bits), min(a[i][i] for i in completions))
                    complement = [states[i] for i in completions if states[i] not in retained]
                    if complement:
                        self.assertLessEqual(oracle.branch_lower(mask, bits, retained), min(rows[s] for s in complement))

    def test_sparse_blocks_leakage_upper_and_complete_coverage(self):
        for name in ('rectangle', 'rectangle_residual', 'square_residual'):
            p = ROOT / 'results/marginal_sparse_molecular' / name
            c = json.loads((p / 'certificate.json').read_text())
            states, a = matrix(c)
            oracle = DeterminantOracle(c)
            retained = c['retained_states']
            indices = [states.index(s) for s in retained]
            expected = prepare(a, indices)
            actual = oracle.retained_data(retained, F(c['complement_lower']))
            self.assertEqual(actual[:2], expected[:2])
            self.assertLessEqual(actual[2], expected[2])
            tree, stats = oracle.cover(retained, actual[2], c['complement_tree'])
            self.assertEqual(tree, c['complement_tree'])
            self.assertEqual(stats['exact_Q_rows'] + stats['pruned_Q_states'], 70 - len(retained))
            witness = c['independent_upper']
            vector = dict(zip(witness['states'], witness['amplitudes']))
            energy = sum(vector.get(s, 0) * vector.get(t, 0) * a[i][j]
                         for i, s in enumerate(states) for j, t in enumerate(states))
            self.assertEqual(oracle.upper(witness), energy / sum(x*x for x in vector.values()))

    def test_saved_replays_never_call_full_matrix_builder(self):
        for name in ('rectangle', 'rectangle_residual', 'square_residual'):
            p = ROOT / 'results/marginal_sparse_molecular' / name
            c = json.loads((p / 'certificate.json').read_text())
            with patch('experiments.marginal_molecular_gershgorin.matrix', side_effect=AssertionError('Full matrix')):
                r = replay(c)
            self.assertEqual(r['width'], json.loads((p / 'receipt.json').read_text())['width'])
            self.assertLess(r['unique_action_states'], 70)
            self.assertGreaterEqual(r['referenced_determinants'], r['unique_action_states'])
            self.assertGreater(r['complement_coverage']['pruned_Q_states'], 0)
            if name.startswith('rectangle'):
                self.assertLess(r['width_float'], 2e-10)
            else:
                self.assertGreater(r['width_float'], 1e-7)

    def test_construction_needs_only_hamiltonian_and_small_retained_eigensolves(self):
        import numpy as np
        from experiments.marginal_sparse_molecular import run
        source = json.loads((ROOT / SOURCES['rectangle']).read_text())
        bare = {key: source[key] for key in ('modes', 'particles', 'hamiltonian')}
        original = np.linalg.eigh
        dimensions = []
        def small_eigh(a):
            dimensions.append(len(a))
            self.assertLessEqual(len(a), 32)
            return original(a)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'input.json').write_text(json.dumps(bare))
            with patch('experiments.marginal_sparse_molecular.ROOT', root), \
                 patch('experiments.marginal_sparse_molecular.SOURCES', {'rectangle': 'input.json'}), \
                 patch('experiments.marginal_molecular_gershgorin.matrix', side_effect=AssertionError('Full matrix')), \
                 patch('numpy.linalg.eigh', side_effect=small_eigh), contextlib.redirect_stdout(io.StringIO()):
                run('rectangle')
            result = root / 'results/marginal_sparse_molecular/rectangle'
            receipt = json.loads((result / 'receipt.json').read_text())
            self.assertLess(receipt['width_float'], 2e-10)
            self.assertLess(receipt['construction_unique_action_states'], 70)
            self.assertEqual(dimensions[0], 1)

    def test_false_coverage_gap_and_witness_are_rejected(self):
        p = ROOT / 'results/marginal_sparse_molecular/rectangle/certificate.json'
        c = json.loads(p.read_text())
        for update in ({'complement_tree': ['bound']}, {'complement_tree': ['split', ['bound']]},
                       {'complement_tree': None}, {'complement_lower': '10'}, {'lower': '10'},
                       {'retained_states': [15, 15]}, {'retained_states': [True]},
                       {'independent_upper': {'states': [15], 'amplitudes': [0]}}):
            bad = copy.deepcopy(c)
            bad.update(update)
            with self.assertRaises(ValueError):
                replay(bad)
        bad = copy.deepcopy(c)
        bad['hamiltonian'] = encode(add(DeterminantOracle(c).h, mono((), F(-1))))
        with self.assertRaises(ValueError): replay(bad)

    def test_exponential_diagonal_sector_covered_without_any_Q_action(self):
        modes, particles = 64, 32
        h = add(*(mono(((1, i), (0, i)), F(i)) for i in range(modes)))
        oracle = DeterminantOracle({'modes': modes, 'particles': particles, 'hamiltonian': encode(h)})
        retained = [(1 << particles) - 1]
        gamma = F(particles * (particles - 1), 2) + F(1, 2)
        with patch.object(oracle, 'action', side_effect=AssertionError('No individual Q action needed')):
            tree, stats = oracle.cover(retained, gamma)
            self.assertEqual(oracle.cover(retained, gamma, tree)[1], stats)
        self.assertEqual(stats['nodes'], 2 * particles + 1)
        self.assertEqual(stats['exact_Q_rows'], 0)
        self.assertEqual(stats['pruned_Q_states'], comb(modes, particles) - 1)


if __name__ == '__main__':
    unittest.main()
