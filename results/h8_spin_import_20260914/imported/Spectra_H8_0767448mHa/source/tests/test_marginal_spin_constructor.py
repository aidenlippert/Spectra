import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_spin_constructor import build, spin_states
from experiments.marginal_spin_reduction import SpinZeroOracle, replay
from experiments.marginal_symbolic import add, mono, encode


def hopping_fixture():
    h = add(*(mono(((1, i+spin), (0, j+spin)), value)
              for spin in (0, 1) for i, j, value in ((0, 2, 1), (2, 0, 1), (2, 2, 2))))
    return {'modes': 4, 'particles': 2, 'hamiltonian': encode(h)}


class DirectSpinConstructorTests(unittest.TestCase):
    def test_saved_bare_h6_interval_and_provenance(self):
        root = Path(__file__).resolve().parents[1]
        path = root / 'results/marginal_h6/direct_spin'
        c = json.loads((path / 'certificate.json').read_text())
        original = json.loads((root / 'results/marginal_h6/fixture.json').read_text())
        for key in ('modes', 'particles', 'hamiltonian'):
            self.assertEqual(c[key], original[key])
        r = replay(c)
        saved = json.loads((path / 'receipt.json').read_text())
        self.assertEqual(r['width'], saved['width'])
        self.assertLess(r['width_float'], 2e-10)
        self.assertEqual(r['response_dimension'], 28)
        self.assertEqual(r['unique_determinant_sources'], 400)
        self.assertEqual(saved['imported_fields'], ['modes', 'particles', 'hamiltonian'])
        self.assertEqual(saved['construction_unique_determinants'], 400)
        selection = json.loads((path / 'selection.json').read_text())
        self.assertEqual(c['spin_symmetric_certificate']['retained_states'], selection['retained_states'])

    def test_bare_constructor_ignores_poisoned_proofs_and_uses_only_reduced_states(self):
        import numpy as np
        from itertools import combinations
        source_action, read_text = DeterminantOracle.action, Path.read_text
        eigh, eigvalsh = np.linalg.eigh, np.linalg.eigvalsh
        sources = set()
        def action(oracle, state):
            self.assertEqual(state.bit_count(), 2)
            self.assertEqual((state & 5).bit_count(), 1)
            sources.add(state)
            return source_action(oracle, state)
        def small_eigh(a):
            self.assertLessEqual(len(a), 3)
            return eigh(a)
        def small_eigvalsh(a):
            self.assertLessEqual(len(a), 3)
            return eigvalsh(a)
        def reduced_combinations(items, n):
            self.assertEqual(list(items), [0, 1])
            self.assertEqual(n, 1)
            return combinations(items, n)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            clean = root / 'bare.json'
            poisoned = root / 'poisoned.json'
            base = hopping_fixture()
            clean.write_text(json.dumps(base))
            poisoned.write_text(json.dumps(dict(base, retained_states=[5], complement_lower='100000',
                independent_upper={'states': [5], 'amplitudes': [1]}, blocks=[],
                response_basis='poison', spin_symmetric_certificate='poison')))
            certificates = []
            for source, name in ((clean, 'clean_out'), (poisoned, 'poisoned_out')):
                out = root / name
                def bounded_read(path, *args, **kwargs):
                    self.assertTrue(path == source or out in path.parents, str(path))
                    return read_text(path, *args, **kwargs)
                with patch.object(DeterminantOracle, 'action', action), patch.object(Path, 'read_text', bounded_read), \
                     patch('numpy.linalg.eigh', side_effect=small_eigh), patch('numpy.linalg.eigvalsh', side_effect=small_eigvalsh), \
                     patch('experiments.marginal_spin_constructor.combinations', side_effect=reduced_combinations), \
                     patch('experiments.marginal_spin_reduction.reduce_certificate', side_effect=AssertionError('Inherited proof')), \
                     contextlib.redirect_stdout(io.StringIO()):
                    build(source, out, reference_cap=1, upper_steps=3)
                certificates.append(json.loads((out / 'certificate.json').read_text()))
                r = json.loads((out / 'receipt.json').read_text())
                self.assertLess(r['width_float'], 1e-7)
                self.assertEqual(r['imported_fields'], ['modes', 'particles', 'hamiltonian'])
                self.assertEqual(r['construction_unique_determinants'], 4)
                self.assertEqual(r['construction_unique_source_determinants'], 4)
                self.assertEqual(r['largest_numerical_eigensolve'], 3)
                self.assertEqual(replay(certificates[-1])['width'], r['width'])
            self.assertEqual(certificates[0], certificates[1])
            self.assertEqual(sources, {3, 6, 9, 12})
            with self.assertRaises(ValueError): build(clean, root / 'clean_out')

    def test_disconnected_bad_reference_is_not_accepted(self):
        h = add(*(mono(((1, i), (0, i)), 2 if i < 2 else -1) for i in range(4)))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.json'
            source.write_text(json.dumps({'modes': 4, 'particles': 2, 'hamiltonian': encode(h)}))
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
                build(source, root / 'out', reference_cap=1, upper_steps=3)
            self.assertFalse((root / 'out/certificate.json').exists())
            self.assertEqual(json.loads((root / 'out/failure.json').read_text())['status'], 'not_accepted')

    def test_explicit_sector_budget_rejects_before_enumeration(self):
        oracle = SpinZeroOracle({'modes': 16, 'particles': 8, 'hamiltonian': []})
        with self.assertRaises(ValueError): spin_states(oracle)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.json'
            source.write_text(json.dumps({'modes': 16, 'particles': 8, 'hamiltonian': []}))
            with self.assertRaises(ValueError): build(source, root / 'out')
            self.assertFalse((root / 'out').exists())


if __name__ == '__main__':
    unittest.main()
