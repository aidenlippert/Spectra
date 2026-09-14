import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_sparse_upper import replay, refine
from experiments.marginal_symbolic import add, encode, mono

ROOT = Path(__file__).resolve().parents[1]


class SparseUpperTests(unittest.TestCase):
    def test_exact_rayleigh_residual_and_witness_validation(self):
        h = add(mono(((1, 0), (0, 1))), mono(((1, 1), (0, 0))), mono(((1, 1), (0, 1)), F(2)))
        c = {'kind': 'sparse_variational_upper_v1', 'modes': 2, 'particles': 1,
             'hamiltonian': encode(h), 'independent_upper': {'states': [1, 2], 'amplitudes': [1, 1]}}
        r = replay(c)
        self.assertEqual(F(r['upper']), F(2))
        self.assertEqual(F(r['residual_norm_squared']), F(1))
        bad = copy.deepcopy(c)
        bad['independent_upper']['amplitudes'] = [0, 0]
        with self.assertRaises(ValueError): replay(bad)

    def test_saved_upper_replay_and_monotone_history(self):
        for steps in (12, 20):
            p = ROOT / f'results/marginal_h6/krylov_upper_{steps}'
            c = json.loads((p / 'certificate.json').read_text())
            with patch('experiments.marginal_molecular_gershgorin.matrix', side_effect=AssertionError('Full sector matrix')):
                r = replay(c)
            expected = json.loads((p / 'receipt.json').read_text())
            self.assertEqual(r['upper'], expected['upper'])
            self.assertEqual(r['residual_norm_squared'], expected['residual_norm_squared'])
            self.assertEqual(r['unique_action_states'], 200)
            history = json.loads((p / 'history.json').read_text())
            energies = [F(x['upper']) for x in history]
            self.assertTrue(all(b <= a for a, b in zip(energies, energies[1:])))
            self.assertLess(F(r['upper']), F(expected['initial_upper']))

    def test_physical_krylov_proposer_on_small_full_response(self):
        import numpy as np
        h = add(mono(((1, 0), (0, 1))), mono(((1, 1), (0, 0))), mono(((1, 1), (0, 1)), F(2)))
        fixture = {'modes': 2, 'particles': 1, 'hamiltonian': encode(h)}
        proposal = {'independent_upper': {'states': [1], 'amplitudes': [1]}}
        original = np.linalg.eigh
        def small(a):
            self.assertLessEqual(len(a), 2)
            return original(a)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'fixture.json').write_text(json.dumps(fixture))
            (root / 'proposal.json').write_text(json.dumps(proposal))
            with patch('numpy.linalg.eigh', side_effect=small), contextlib.redirect_stdout(io.StringIO()):
                refine(root / 'fixture.json', root / 'proposal.json', root / 'out', 2)
            r = json.loads((root / 'out/receipt.json').read_text())
            self.assertLess(r['upper_float'], -.4142)
            self.assertLess(r['residual_norm_squared_float'], 1e-18)


if __name__ == '__main__':
    unittest.main()
