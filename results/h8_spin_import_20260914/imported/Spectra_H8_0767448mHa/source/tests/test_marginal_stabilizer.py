import unittest
import numpy as np

from experiments.marginal_stabilizer import split_blocks


class StabilizerTests(unittest.TestCase):
    def test_selected_dictionary_split_is_complete_and_validated(self):
        blocks,receipt=split_blocks(4,4,indices=[0])
        self.assertEqual(len(receipt),1)
        self.assertTrue(all(b['name'].split(':')[0]=='0' for b in blocks))
        u=np.column_stack([b['basis_transform'] for b in blocks])
        self.assertLess(np.max(abs(u@u.T-np.eye(u.shape[0]))),1e-12)
        with self.assertRaises(ValueError):split_blocks(4,4,indices=[-1])
        with self.assertRaises(ValueError):split_blocks(4,4,indices=[False])

    def test_subspaces_cover_each_original_dictionary(self):
        blocks,receipt=split_blocks(4,4)
        grouped={}
        for b in blocks:grouped.setdefault(b['name'].split(':')[0],[]).append(b['basis_transform'])
        for index,r in enumerate(receipt):
            u=np.column_stack(grouped[str(index)])
            self.assertEqual(u.shape,(r['original_dimension'],r['original_dimension']))
            self.assertLess(np.max(np.abs(u@u.T-np.eye(len(u)))),1e-12)
            self.assertEqual(sum(r['subspaces']),r['original_dimension'])


if __name__=='__main__':unittest.main()
