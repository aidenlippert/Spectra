from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import add, adj, canonical, mono, product
from research.interacting_scaling_20260915.dictionary import clustered_frame, close_rows, ideal_basis, representative, windows
from research.sector_quotient_20260914.fast_twirl import twirl


class DirectDictionaryTests(unittest.TestCase):
    def test_nested_blocks_keep_long_range_cross_terms(self):
        groups, specs = clustered_frame(8, windows(4, 2))
        extended, ext_specs = clustered_frame(8, windows(4, 2)+windows(4, 3))
        self.assertEqual(groups, extended[:len(groups)])
        for old, new in zip(specs, ext_specs):
            self.assertEqual(old[0], new[0])
            self.assertTrue((old[1] == new[1]).all())
        collective = [g for g in groups if g['name'].startswith('collective_pair_supports:mixed-')]
        self.assertTrue(any(any({i//2 for c, i in w} == {0, 3} for w in g['words']) for g in collective))

    def test_rows_retain_exact_hermitian_products_and_spin_average(self):
        left = add(mono(((1, 0), (0, 3), (0, 2)), 2), mono(((0, 4),), -3))
        square = product(adj(left), left)
        rows = set(close_rows(square))
        self.assertTrue(all(representative(w) in rows for w in square))
        self.assertTrue(all(representative(w) in rows for w in twirl(square)))
        for w in rows:
            p = canonical(mono(w))
            self.assertTrue(all(representative(v) in rows for v in twirl(add(p, canonical(adj(p))))))

    def test_ideals_allow_global_charge_transfer(self):
        basis = ideal_basis(8, windows(4, 2))
        self.assertIn(canonical(mono(((1, 0), (1, 7), (0, 7), (0, 0)))), basis)
        transfer = {( (1, 0), (0, 6)): F(1, 2), ((1, 6), (0, 0)): F(1, 2)}
        self.assertIn(transfer, basis)

    def test_invalid_clusters_and_nonconserving_rows_refused(self):
        for clusters in [[(0, 0)], [(0, 5)], [(0,)]]:
            with self.assertRaises(ValueError):
                clustered_frame(8, clusters)
        with self.assertRaises(ValueError):
            representative(((1, 0),))

    def test_collective_cluster_extension_contains_old_spin_space(self):
        import numpy as np
        from research.interacting_scaling_20260915.warmstart import basis_embedding
        groups, specs = clustered_frame(6, [(0, 1), (1, 2)])
        larger, new_specs = clustered_frame(6, [(0, 1), (1, 2)], collective_clusters=[(0, 1, 2)])
        by_name = {g['name']: g for g in larger}
        self.assertTrue({g['name'] for g in groups} <= set(by_name))
        for old in groups:
            new = by_name[old['name']]
            self.assertTrue(set(old['words']) <= set(new['words']))
        new_by_name = {larger[g]['name']: (g, W) for g, W, _ in new_specs}
        for g, V, _ in specs:
            new_g, W = new_by_name[groups[g]['name']]
            embedding = np.zeros((len(larger[new_g]['words']), V.shape[1]))
            positions = {w: i for i, w in enumerate(larger[new_g]['words'])}
            for i, word in enumerate(groups[g]['words']): embedding[positions[word]] = V[i]
            recovered = W@np.linalg.lstsq(W, embedding, rcond=None)[0]
            self.assertLess(np.max(abs(recovered-embedding)), 1e-11)
            A = basis_embedding(groups[g]['words'], V, larger[new_g]['words'], W)
            self.assertLess(np.max(abs(W@A-embedding)), 1e-14)
        with self.assertRaises(ValueError):
            basis_embedding([((0, 7),)], np.eye(1), [((0, 0),)], np.eye(1))


if __name__ == '__main__':
    unittest.main()
