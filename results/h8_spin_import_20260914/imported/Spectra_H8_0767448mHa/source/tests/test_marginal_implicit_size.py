from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import copy
import json
import unittest

from experiments.marginal_defect_dicke import DefectDicke, reference_complement, enlarged_workspace, workspace_data
from experiments.marginal_implicit_certificate import certificate_modes, replay
from experiments.marginal_general_schur import fixture, apply_columns, gram
from experiments.marginal_schur_transfer import complement_lower
from experiments.marginal_sector_reference import above_lower_bound, jacobi_data
from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_symbolic import decode, transform
from tests.test_marginal_defect_dicke import expand

ROOT = Path(__file__).resolve().parents[1]


class ImplicitSizeTests(unittest.TestCase):
    def test_complement_bound_covers_every_charge_and_spin_block(self):
        for m in range(2, 9):
            c = reference_complement(m)
            total = 0
            for d in range(m // 2 + 1):
                s = m - 2 * d
                charge_copies = factorial(m) // (factorial(d) ** 2 * factorial(s))
                for k in range(s // 2 + 1):
                    spin_copies = comb(s, k) - (comb(s, k - 1) if k else 0)
                    twice_j = s - 2 * k
                    diagonal = [F(m * (m - 2), 4) + (F(-twice_j, 2) + q) ** 2 for q in range(twice_j + 1)]
                    squared = [F((q + 1) * (twice_j - q), 25) for q in range(twice_j)]
                    total += charge_copies * spin_copies * (twice_j + 1)
                    if d == 1 and k == 0:
                        self.assertEqual((diagonal, squared), jacobi_data(2 * m, F(1, 5), 1))
                    if d or k:
                        self.assertTrue(above_lower_bound(diagonal, squared, c), (m, d, k))
            self.assertEqual(total, comb(2 * m, m))
        self.assertEqual(reference_complement(5), complement_lower())

    def test_generalized_fixtures_preserve_old_model_and_break_new_charges(self):
        old = json.loads((ROOT / 'results/marginal_enlarged_schur/mixed_1_100/certificate.json').read_text())
        self.assertEqual(fixture(F(1, 100), True, 10), decode(old['hamiltonian'], 10, 4))
        for mixed in (False, True):
            h = fixture(F(1, 100), mixed, 12)
            for pair in range(6):
                self.assertTrue(any(sum((1 if c else -1) for c, i in word if i in (pair, pair + 6)) for word in h))
            self.assertNotEqual(h, transform(h, [(i + 6) % 12 for i in range(12)]))
        with self.assertRaises(ValueError): fixture(modes=11)
        with self.assertRaises(ValueError): fixture(interaction=True, modes=4)

    def test_twelve_mode_embedding_and_projected_operators_match_fock_oracle(self):
        h = fixture(F(1, 100), modes=12)
        workspace = enlarged_workspace(h, 6)
        self.assertEqual(len(workspace['basis']), 17)
        model = workspace['model']
        z = model.reference()
        h0 = hopping_polynomial(12, F(1, 5))
        self.assertEqual([expand(model, model.reference_action(v)) for v in z], apply_columns(h0, [expand(model, v) for v in z]))
        u = [expand(model, v) for v in workspace['basis']]
        hu = apply_columns(h, u)
        self.assertEqual(hu, [expand(model, v) for v in workspace['action']])
        data = workspace_data(workspace)
        self.assertEqual(data['metric'], gram(u, u))
        self.assertEqual(data['projected_h'], gram(u, hu))
        self.assertTrue(all(x == 0 for row in data['leakage'][:7] for x in row))

    def test_twelve_mode_saved_energy_interval(self):
        certificate = json.loads((ROOT / 'results/marginal_implicit_certificate/m12_cycle_1_100/certificate.json').read_text())
        receipt = replay(certificate)
        self.assertEqual(receipt['modes'], 12)
        self.assertEqual(receipt['particles'], 6)
        self.assertEqual(receipt['stage_dimensions'], [7, 17])
        self.assertLess(receipt['width_float'], 1e-3)
        bad = copy.deepcopy(certificate)
        bad['lower'] = '10'
        with self.assertRaises(ValueError): replay(bad)

    def test_size_metadata_and_reference_limits(self):
        for modes, particles in ((True, 1), (3, 1), (11, 5), (12, 5), (12, True)):
            with self.assertRaises(ValueError): certificate_modes({'kind': 'fully_implicit_schur_v1', 'modes': modes, 'particles': particles})
        with self.assertRaises(ValueError): reference_complement(1)
        with self.assertRaises(ValueError): enlarged_workspace({}, 64)

    def test_legacy_target_and_explicit_chain_have_identical_replay(self):
        path = ROOT / 'results/marginal_implicit_certificate/cycle_1_100_targeted_compressed/certificate.json'
        certificate = json.loads(path.read_text())
        original = replay(certificate)
        chain = copy.deepcopy(certificate)
        chain['target_chain'] = [chain.pop('target_coefficients')]
        self.assertEqual(replay(chain), original)
        chain['target_coefficients'] = certificate['target_coefficients']
        with self.assertRaises(ValueError): replay(chain)
        chain['target_coefficients'] = None
        chain['target_chain'] = [[], [], [], []]
        with self.assertRaises(ValueError): replay(chain)


if __name__ == '__main__':
    unittest.main()
