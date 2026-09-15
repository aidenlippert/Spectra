from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import unittest
import numpy as np
from experiments.v4_physics import Model, Experiment, make_family, expectation, likelihood_table, sample


def matrix(word):
    matrices = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]),
                'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}
    out = np.eye(1)
    for char in word: out = np.kron(out, matrices[char])
    return out


class PhysicsTests(unittest.TestCase):
    def test_all_hypotheses_and_ordered_gates_against_dense_born_rule(self):
        for n in (2, 3, 4):
            for prefix in (False, True):
                models, actions = make_family(n, seed=2, prefix=prefix)
                self.assertEqual({m.target_sign for m in models}, {-1, 1})
                self.assertLessEqual(len(actions), 12)
                for model in models:
                    unitary = np.eye(2 ** n, dtype=complex)
                    for p, sign in model.gates:
                        unitary = (np.eye(2 ** n) - 1j * sign * matrix(p)) / np.sqrt(2) @ unitary
                    for action in actions:
                        rho = (np.eye(2 ** n) + matrix(action.prep)) / 2 ** n
                        dense = np.trace(rho @ unitary.conj().T @ matrix(action.measurement) @ unitary).real
                        self.assertAlmostEqual(dense, expectation(model, action), places=12)
                        self.assertEqual(expectation(model, action), expectation(replace(model, target_sign=-model.target_sign), action))

    def test_exact_channel_and_each_candidate_has_sign_information(self):
        models, actions = make_family(3, 31, True)
        table = likelihood_table(models, actions, (F(3, 5),) * len(actions))
        self.assertTrue(all(type(p) is F and 0 <= p <= 1 for row in table for p in row))
        for i in range(0, 8, 2): self.assertNotEqual(table[i], table[i + 1])
        rng = np.random.default_rng(1)
        self.assertTrue(set(sample(models[0], actions[0], 20, F(3, 5), rng)) <= {0, 1})

    def test_invalid_physical_inputs(self):
        with self.assertRaises(ValueError): Model((('XX', True),), 1)
        with self.assertRaises(ValueError): Model((('XYZ', 1),), 1)
        with self.assertRaises(ValueError): Experiment('II', 'XI', 1)
        with self.assertRaises(ValueError): make_family(5)
        models, actions = make_family(2)
        with self.assertRaises(ValueError): likelihood_table(models, actions, (0.5,) * len(actions))


if __name__ == '__main__': unittest.main()
