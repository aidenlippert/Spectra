"""The whole-frame ceiling must cover every discovery generator and new gate."""
from pathlib import Path
import json
import unittest

from research.spin_enrichment_20260913.full_dual import check, full_groups


class FullFamily(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from research.spin_completion_20260913.discovery import Model
        root = Path(__file__).resolve().parents[2]
        prior = root/'results/molecular_collective_20260913/campaign/h6'
        cls.model = Model(json.loads((prior/'fixture.json').read_text()), json.loads((prior/'rank_10/tail.json').read_text()))
        cls.previous_dual = json.loads((root/'results/spin_subspace_20260913/full_dual/witness.json').read_text())

    def test_rebuilt_exact_frame_matches_every_numerical_generator(self):
        m = self.model; base, added = full_groups(m.p, m.tail, m.symmetry['parity_masks'])
        self.assertEqual([g['polynomials'] for g in base], [g['polynomials'] for g in m.base])
        self.assertEqual([g['polynomials'] for g in added], [g['polynomials'] for g in m.frames])
        self.assertEqual(sum(len(g['polynomials']) for g in added), 492)

    def test_previous_diagonal_spin_dual_fails_the_added_PSD_gate(self):
        witness = {**self.previous_dual, 'kind': 'spin_completion_full_dual_v1'}
        with self.assertRaisesRegex(ValueError, 'Negative exact PSD pivot'):
            check(self.model.data, self.model.tail, witness)


if __name__ == '__main__':
    unittest.main()
