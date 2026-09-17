import unittest
import numpy as np
from research.constructive_response_20260916.interacting_probe import (
    balanced_basis,local_rotation,lifted_gate,
)
from research.constructive_response_20260916.layered_probe import fermion_gate,permutation_sign


class FermionEmbeddingTests(unittest.TestCase):
    def test_adjacent_embedding_matches_original(self):
        labels=balanced_basis(4);g=local_rotation()
        a=fermion_gate(labels,4,2,3,g).toarray()
        b=lifted_gate(labels,1,g).toarray()
        np.testing.assert_array_equal(a,b)

    def test_nonadjacent_gates_are_orthogonal_and_even_local_gates_commute(self):
        labels=balanced_basis(4);g=local_rotation()
        a=fermion_gate(labels,4,0,2,g).toarray()
        b=fermion_gate(labels,4,1,3,g).toarray()
        np.testing.assert_allclose(a.T@a,np.eye(len(labels)),atol=5e-16)
        np.testing.assert_allclose(a@b,b@a,atol=5e-16)

    def test_occupied_permutation_has_fermion_sign(self):
        self.assertEqual(permutation_sign(0b11,[1,0,2,3]),-1)
        self.assertEqual(permutation_sign(0b01,[1,0,2,3]),1)


if __name__=='__main__':
    unittest.main()
