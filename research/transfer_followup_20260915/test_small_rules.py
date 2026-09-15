from copy import deepcopy
from fractions import Fraction as F
import unittest
from research.transfer_solver_20260915.independent_oracle import matrix, positive_shift
from research.transfer_followup_20260915.enumerated_baseline import blocks, quotient, shifted
from research.transfer_followup_20260915.fragments import validate_decomposition


class SmallRules(unittest.TestCase):
    def test_independent_fermion_sign(self):
        data = {'modes': 2, 'particles': 2, 'hamiltonian': [
            {'word': [[1, 1], [1, 0], [0, 1], [0, 0]], 'coefficient': '1'}]}
        states, H = matrix(data)
        self.assertEqual(states, [3])
        self.assertEqual(H, [[F(-1)]])
        self.assertTrue(positive_shift(H, F(-2))[0])
        self.assertFalse(positive_shift(H, F(0))[0])

    def test_all_number_blocks_and_exact_quotient(self):
        data = {'modes': 4, 'particles': 2, 'hamiltonian': [
            {'word': [[1, i], [0, i]], 'coefficient': str(i+1)} for i in range(4)]}
        groups, den, cost = blocks(data)
        self.assertEqual(cost['full_sector_labels_enumerated'], 6)
        self.assertEqual(sum(len(states) for a, states, H in groups), 6)
        for alpha, states, H in groups:
            for j, determinant in enumerate(states):
                vector = [int(i == j) for i in range(len(states))]
                expected = sum(i+1 for i in range(4) if determinant & (1 << i))
                self.assertEqual(quotient(H, den, vector), expected)
                K, kd = shifted(H, den, F(1, 3))
                self.assertEqual(F(K[j][j], kd), F(expected)-F(1, 3))

    def test_fragment_sector_and_interaction_are_required(self):
        fragment = {'modes': 2, 'particles': 1, 'hamiltonian': [
            {'word': [[1, 0], [0, 0]], 'coefficient': '1'}]}
        data = {'modes': 4, 'particles': 2, 'hamiltonian': fragment['hamiltonian']+[
            {'word': [[1, 2], [0, 2]], 'coefficient': '1'}],
            'fragment_particle_constraints': [{'orbitals': [0, 1], 'particles': 1}, {'orbitals': [2, 3], 'particles': 1}]}
        validate_decomposition(data, fragment, 2)
        missing = deepcopy(data)
        missing.pop('fragment_particle_constraints')
        with self.assertRaisesRegex(ValueError, 'Local particle'):
            validate_decomposition(missing, fragment, 2)
        coupled = deepcopy(data)
        coupled['hamiltonian'].append({'word': [[1, 0], [0, 2]], 'coefficient': '1'})
        with self.assertRaisesRegex(ValueError, 'disjoint sum'):
            validate_decomposition(coupled, fragment, 2)


if __name__ == '__main__':
    unittest.main()
