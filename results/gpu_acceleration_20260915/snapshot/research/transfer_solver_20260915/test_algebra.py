"""Independent tiny-algebra checks for input-derived spin grouping."""
import unittest
from fractions import Fraction as F
import numpy as np
from experiments.marginal_symbolic import mono, add, product, scale, canonical
from research.transfer_solver_20260915.algebra import frame, raising, highest_weights, normalize


class SpinGroups(unittest.TestCase):
    def test_all_mixed_words_against_literal_car_and_reordered_groups(self):
        for m in (4, 6, 8, 12):
            h = add(*(mono(((1, i), (0, i)), F(i//2+1)) for i in range(m)))
            groups, _, _, _ = frame(h, m)
            groups = list(reversed(groups))
            ids, matrices, targets = raising(groups)
            plus = add(*(mono(((1, i), (0, i+1))) for i in range(0, m, 2)))
            for g in ids:
                for j, word in enumerate(groups[g]['words']):
                    expected = add(product(plus, mono(word)), scale(product(mono(word), plus), -1))
                    got = {}
                    if g in matrices:
                        col = matrices[g].getcol(j).tocoo()
                        got = canonical(add(*(mono(groups[targets[g]]['words'][i], F(int(v)))
                                              for i, v in zip(col.row, col.data))))
                    self.assertEqual(got, expected)
            specs = highest_weights(groups)
            self.assertTrue(all(v.shape[1] for _, v, _ in specs))

    def test_asymmetric_hopping_removes_spatial_parity_without_fixed_indices(self):
        m = 8
        diag = add(*(mono(((1, i), (0, i)), F(i//2+1)) for i in range(m)))
        hop = add(*(mono(((1, i), (0, i+2))) for i in range(m-2)),
                  *(mono(((1, i+2), (0, i))) for i in range(m-2)))
        g0, _, _, s0 = frame(diag, m)
        g1, _, _, s1 = frame(add(diag, hop), m)
        self.assertGreater(len(s0['spatial_parity_masks']), len(s1['spatial_parity_masks']))
        self.assertNotEqual(len(g0), len(g1))
        self.assertTrue(highest_weights(g1))

    def test_missing_spin_partner_refused(self):
        m = 6
        h = add(*(mono(((1, i), (0, i))) for i in range(m)))
        groups, _, _, _ = frame(h, m)
        ids, matrices, targets = raising(groups)
        target = next(iter(targets.values()))
        groups[target]['words'] = groups[target]['words'][1:]
        with self.assertRaises((ValueError, IndexError)):
            highest_weights(groups)


if __name__ == '__main__':
    unittest.main()
