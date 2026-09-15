from fractions import Fraction as F
from itertools import product
import unittest
from unittest.mock import patch
import numpy as np
from experiments.v4_calibration import (learn_from_calibration, labelled_source, test_world,
    sample_budget, exact_policy_risk, policy_from_estimated_moments)
from experiments.v4_moments import signed_gram_plan


class CalibrationTests(unittest.TestCase):
    def test_data_only_learning_and_costed_confidence(self):
        p, L, s = test_world(True)
        source = labelled_source(L, s, 22001)
        with patch('experiments.v4_calibration.signed_gram_plan', side_effect=AssertionError('hidden teacher access')):
            learned = learn_from_calibration(source, 4)
        actual = exact_policy_risk(learned, p, L, s)
        self.assertEqual(actual, F(9, 50))
        lo, hi = learned['conditional_risk_interval']
        self.assertLessEqual(lo, actual); self.assertLessEqual(actual, hi)
        self.assertLessEqual(learned['union_failure_upper_bound'], learned['delta'])
        self.assertEqual(learned['work']['labelled_episodes'], 15 * 38400)
        self.assertEqual(learned['work']['visible_readouts'], 24 * 38400)
        self.assertEqual(learned['work']['target_label_readouts'], 15 * 38400)

    def test_repeat_actions_require_same_model_independent_noise(self):
        L, s = ((F(9, 10),), (F(1, 10),)), (1, 1)
        labels, outcomes = labelled_source(L, s, 23001)((0, 0), 30000)
        gram = float(np.mean(labels * outcomes[:, 0] * outcomes[:, 1]))
        self.assertLess(abs(gram - .41), .015)
        labels, outcomes = labelled_source(L, s, 23001, redraw_each_action=True)((0, 0), 30000)
        broken = float(np.mean(labels * outcomes[:, 0] * outcomes[:, 1]))
        self.assertLess(abs(broken - .25), .015)
        self.assertGreater(abs(gram - broken), .12)

    def test_uniform_surrogate_bound_all_small_action_trees(self):
        rng = np.random.default_rng(24001)
        for _ in range(8):
            p = (F(1, 4),) * 4; signs = (1, 1, -1, -1)
            L = tuple(tuple(F(int(x), 5) for x in row) for row in rng.integers(0, 6, size=(4, 2)))
            exact = signed_gram_plan(p, L, signs)
            mu = exact['moments']; epsilon = F(1, 20)
            bias = mu['bias'] + epsilon
            first = tuple(x - epsilon for x in mu['first'])
            gram = tuple(tuple(mu['gram'][a][b] + (epsilon if a == b else -epsilon) for b in range(2)) for a in range(2))
            learned = policy_from_estimated_moments(bias, first, gram, epsilon)
            kappa = F(9, 2) * epsilon
            for a, b0, b1 in product(range(2), repeat=3):
                for leaves in product((-1, 1), repeat=4):
                    policy = {'first_action': a, 'child_actions': (b0, b1),
                              'leaf_decisions': (leaves[:2], leaves[2:])}
                    actual = exact_policy_risk(policy, p, L, signs)
                    masses = (bias-first[a]-first[b0]+gram[a][b0], first[b0]-gram[a][b0],
                              first[a]-gram[a][b1], gram[a][b1])
                    surrogate = F(1, 2) - sum(d * mass for d, mass in zip(leaves, masses)) / 2
                    self.assertLessEqual(abs(actual-surrogate), kappa)
                    self.assertLessEqual(learned['surrogate_risk'], surrogate)
            self.assertLessEqual(exact_policy_risk(learned,p,L,signs)-min(exact['values']), 2*kappa)

    def test_invalid_records_budget_and_inconsistent_moment_refusal(self):
        with self.assertRaises(ValueError): sample_budget(4, F(1, 10000), F(1, 100))
        with self.assertRaises(ValueError): sample_budget(True, F(1, 40), F(1, 100))
        def malformed(actions,n): return np.ones(n), np.zeros((n, len(actions)+1))
        with self.assertRaises(ValueError): learn_from_calibration(malformed, 1)
        bad = policy_from_estimated_moments(F(0), (F(0),), ((F(1),),), F(1, 100))
        self.assertFalse(bad['accepted'])
        self.assertIsNone(bad['conditional_risk_interval'])


if __name__ == '__main__': unittest.main()
