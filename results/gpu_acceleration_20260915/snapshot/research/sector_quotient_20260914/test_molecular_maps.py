"""Integration checks on the actual H8 full-degree operator maps."""
import json
import unittest
import numpy as np
from scipy import sparse
from research.sector_quotient_20260914.budget import OUT


class MolecularMapTests(unittest.TestCase):
    def test_all_paired_sextic_maps_cancel_with_transposed_indices(self):
        prepared = OUT/'prepared'; meta = json.loads((prepared/'frame.json').read_text())
        selected = np.load(prepared/'selected.npy')
        sixth = np.array([len(meta['rows'][i]) == 6 for i in selected])
        rng = np.random.default_rng(1921)
        for pair in meta['pairs']:
            if len(pair['members']) != 2: continue
            i, j = pair['members']; n = len(meta['groups'][i]['words'])
            M = sparse.load_npz(prepared/f'map_{i}.npz'); B = sparse.load_npz(prepared/f'map_{j}.npz')
            order = np.array(pair['adjoint_order'])
            paired = M+B[:, (order[None, :]*n+order[:, None]).ravel()]
            paired.eliminate_zeros()
            self.assertEqual(paired[sixth, :].nnz, 0)
            q = rng.normal(size=(n, n)); q = (q+q.T)/2
            other = np.empty_like(q); other[np.ix_(order, order)] = q
            self.assertTrue(np.allclose(paired@q.ravel(), M@q.ravel()+B@other.ravel(), atol=1e-11))

    def test_disjoint_spatial_triple_constraint_is_retained(self):
        prepared = OUT/'prepared'; meta = json.loads((prepared/'frame.json').read_text())
        free = sparse.load_npz(prepared/'free_full.npz')
        T = sparse.load_npz(prepared/'twirl.npz'); selected = np.load(prepared/'selected.npy')
        found = False
        for i, word in enumerate(meta['rows']):
            if len(word) != 6: continue
            left = {p//2 for c, p in word if c}; right = {p//2 for c, p in word if not c}
            if left.isdisjoint(right) and T[selected, i].nnz:
                self.assertEqual(free.getrow(i).nnz, 0)
                found = True; break
        self.assertTrue(found, 'A nonzero contraction-free molecular coefficient must remain constrained')


if __name__ == '__main__': unittest.main()
