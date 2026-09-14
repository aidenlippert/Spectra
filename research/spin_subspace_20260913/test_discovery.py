"""The cached floating pricing map is checked against exact rational Grams."""
from fractions import Fraction as F
from pathlib import Path
import json
import unittest
import numpy as np

from research.certificate_scaling.commutator_dual_witness import moment_decode
from research.joint_patterns_20260913.dual import integer_grams
from research.spin_subspace_20260913.discovery import Model, DIRECTION_DENOMINATOR
from research.spin_subspace_20260913.full_dual import frames, check as check_full_spin


class PricingMap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        previous = root/'results/molecular_collective_20260913/campaign/h6'
        cls.model = Model(json.loads((previous/'fixture.json').read_text()), json.loads((previous/'rank_10/tail.json').read_text()))
        witness = json.loads((root/'results/joint_patterns_20260913/dual/witness.json').read_text())
        cls.y = moment_decode(witness['moments'], 12)
        cls.raw = {'rows': cls.model.rows, 'values': [float(cls.y.get(w, F(0))) for w in cls.model.rows]}

    def test_float_pricing_matches_exact_anticommutator_gram(self):
        m = self.model; frame = m.frames[0]
        matrices, denominator = integer_grams(frame['polynomials'], [self.y], 'anticommutator')
        exact = np.array([[float(F(c, denominator)) for c in row] for row in matrices[0]])
        C = m.frame_coefficients[0]; values = np.array(self.raw['values'])
        actual = C.T@(m.price_maps[0]@values).reshape(C.shape[0], C.shape[0])@C
        self.assertTrue(np.allclose(actual, exact, atol=2e-14, rtol=2e-14))

    def test_rounded_negative_candidates_are_new_operator_directions(self):
        candidates, _ = self.model.price(self.raw, [])
        self.assertEqual([len(group) for group in candidates], [2, 2, 2, 2])
        first = [entry for group in candidates for entry in group]
        next_candidates, _ = self.model.price(self.raw, first)
        for gid, entries in enumerate(next_candidates):
            C = self.model.frame_coefficients[gid]
            old = [C@(np.array(e['vector'])/DIRECTION_DENOMINATOR) for e in first if e['group'] == gid]
            for entry in entries:
                column = C@(np.array(entry['vector'])/DIRECTION_DENOMINATOR)
                before = np.linalg.matrix_rank(np.column_stack(old), tol=1e-7)
                old.append(column)
                self.assertGreater(np.linalg.matrix_rank(np.column_stack(old), tol=1e-7), before)

    def test_exact_dual_rebuild_contains_every_numerical_spin_generator(self):
        model = self.model
        baseline, added = frames(model.p, model.tail, model.symmetry['parity_masks'])
        self.assertEqual([g['polynomials'] for g in added], [g['polynomials'] for g in model.frames])
        self.assertEqual(sum(len(g['polynomials']) for g in added), 252)
        self.assertEqual([g['polynomials'] for g in baseline], [g['polynomials'] for g in model.base])

    def test_full_spin_gate_rejects_the_old_feasible_spin_summed_witness(self):
        root = Path(__file__).resolve().parents[2]
        witness = json.loads((root/'results/joint_patterns_20260913/dual/witness.json').read_text())
        witness['kind'] = 'full_spin_frame_dual_v1'
        with self.assertRaisesRegex(ValueError, 'Negative exact PSD pivot'):
            check_full_spin(self.model.data, self.model.tail, witness)


if __name__ == '__main__':
    unittest.main()
