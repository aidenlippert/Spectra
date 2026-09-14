import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_fixed_point_ldl import propose
from experiments.marginal_h6_complement import retarget, replay_factor
from experiments.marginal_spin_reduction import SpinZeroOracle, spin_operators, commutator, symmetrize, replay
from experiments.marginal_symbolic import add, mono, product, scale, encode

ROOT = Path(__file__).resolve().parents[1]


def ferromagnet():
    plus, minus, z = spin_operators(4)
    h = scale(add(product(z, z), scale(add(product(plus, minus), product(minus, plus)), F(1, 2))), -1)
    base = {'modes': 4, 'particles': 2, 'hamiltonian': encode(h)}
    oracle = SpinZeroOracle(base)
    witnesses = [{'states': [6, 9], 'amplitudes': [1, s]} for s in (1, -1)]
    witness = min(witnesses, key=oracle.upper)
    assert oracle.upper(witness) == -2
    reduced = dict(base, retained_states=[3, 6, 9], complement_lower='-1', lower='-201/100',
                   blocks=[{'states': [12], 'factor': propose([[1]])}],
                   response_basis=[{'states': [12], 'amplitudes': [1]}], independent_upper=witness)
    return dict(base, kind='spin_reduced_perturbation_interval_v1', spin_symmetric_certificate=reduced)


class SpinReductionTests(unittest.TestCase):
    def test_exact_spin_projection_and_approximate_symmetry_rejection(self):
        h = mono(((1, 0), (0, 0)))
        projected = symmetrize(h, 4)
        self.assertEqual(projected, add(mono(((1, 0), (0, 0)), F(1, 2)), mono(((1, 1), (0, 1)), F(1, 2))))
        plus, _, z = spin_operators(4)
        self.assertFalse(commutator(projected, plus))
        self.assertFalse(commutator(projected, z))
        with self.assertRaises(ValueError): SpinZeroOracle({'modes': 4, 'particles': 2, 'hamiltonian': encode(h)})

    def test_high_spin_ground_state_without_spin_adapted_reference(self):
        import numpy as np
        c = ferromagnet()
        original, reduced = DeterminantOracle(c), SpinZeroOracle(c)
        eigenvalues = []
        for oracle in (original, reduced):
            states = [s for s in range(16) if oracle.valid_state(s)]
            a = np.array([[float(oracle.action(t).get(s, 0)) for t in states] for s in states])
            eigenvalues.append(np.linalg.eigvalsh(a)[0])
        self.assertEqual((original.sector_dimension, reduced.sector_dimension), (6, 4))
        self.assertAlmostEqual(eigenvalues[0], -2)
        self.assertAlmostEqual(eigenvalues[1], -2)
        r = replay(c)
        self.assertEqual(F(r['width']), F(1, 100))
        self.assertEqual(r['unique_determinant_sources'], 4)
        self.assertFalse(reduced.valid_state(5))
        with self.assertRaises(ValueError): reduced.cover([3], F(-1))

    def test_original_spin_breaking_is_transferred_by_recomputed_norm_bound(self):
        c = ferromagnet()
        _, _, z = spin_operators(4)
        h = DeterminantOracle(c).h
        c['hamiltonian'] = encode(add(h, z))
        c['spin_symmetry_error_bound'] = '0'  # This supplied claim must never be trusted.
        r = replay(c)
        self.assertEqual(F(r['spin_symmetry_error_bound']), 2)
        self.assertEqual(F(r['lower']), F(-401, 100))
        self.assertEqual(F(r['upper']), -2)
        # The true full-sector ground energy has moved to -3 outside Sz=0.
        oracle = DeterminantOracle(c)
        self.assertEqual(oracle.action(10)[10], -3)
        self.assertLess(F(r['lower']), -3)

    def test_sector_symmetry_and_coverage_tampering(self):
        c = ferromagnet()
        for update in ({'modes': 5}, {'particles': 1}, {'hamiltonian': encode(mono(((1, 0), (0, 0))))},
                       {'blocks': []}, {'retained_states': [3, 5, 9]}, {'lower': '0'}):
            bad = copy.deepcopy(c)
            bad['spin_symmetric_certificate'].update(update)
            with self.assertRaises(ValueError): replay(bad)
        oracle = SpinZeroOracle(c)
        with self.assertRaises(ValueError): oracle.retained([3, 6, 9, 12])

    def test_sharper_gap_retarget_rechecks_and_preserves_outputs(self):
        c = {'kind': 'factor_complement_bound_v1', 'modes': 2, 'particles': 1,
             'hamiltonian': encode(mono(((1, 1), (0, 1)), 2)),
             'retained_states': [1], 'complement_lower': '0',
             'blocks': [{'states': [2], 'factor': propose([[2]])}]}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'source.json'
            source.write_text(json.dumps(c))
            with contextlib.redirect_stdout(io.StringIO()): retarget(source, root / 'out', '3/2')
            got = json.loads((root / 'out/certificate.json').read_text())
            self.assertEqual(replay_factor(got)['complement_lower'], '3/2')
            with self.assertRaises(ValueError): retarget(source, root / 'out', '1')
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
                retarget(source, root / 'false', '3')
            self.assertFalse((root / 'false').exists())

    def test_saved_h6_reduced_global_interval(self):
        for name, directions in (('spin_reduced', 32), ('spin_reduced_sharp', 28)):
            with self.subTest(name=name):
                path = ROOT / 'results/marginal_h6' / name
                c = json.loads((path / 'certificate.json').read_text())
                r = replay(c)
                saved = json.loads((path / 'receipt.json').read_text())
                self.assertEqual(r['width'], saved['width'])
                self.assertLess(r['width_float'], 2e-10)
                self.assertEqual(F(r['spin_symmetry_error_bound']), F(27, 500000000000))
                self.assertEqual(r['unique_determinant_sources'], 400)
                self.assertEqual((r['spin_symmetric_action_states'], r['original_action_states']), (400, 200))
                self.assertEqual(r['complement_coverage']['block_dimensions'], [200, 168])
                self.assertEqual(r['full_sector_dimension'], 924)
                self.assertEqual(r['response_dimension'], directions)


if __name__ == '__main__':
    unittest.main()
