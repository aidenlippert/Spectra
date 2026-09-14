import unittest
from fractions import Fraction as F
from experiments.v5_bayes import (prior_empty, extend_belief, update_belief,
    choose_probe, current_target_decision, nominal_likelihood)
from experiments.v5_frame import null_probe


class BayesTests(unittest.TestCase):
    def test_magnitude_outcome_identifies_first_mode_both_ways(self):
        prior = extend_belief(prior_empty()); probe = null_probe(())
        for outcome, expected in ((0, 1), (1, -1)):
            post = update_belief(prior, 1, F(1), probe, outcome)
            self.assertEqual(current_target_decision(post), expected)
            self.assertGreater(post[(expected,)], .999)
            self.assertAlmostEqual(sum(post.values()), 1)
        self.assertGreater(nominal_likelihood((1,), F(1), probe, 0), .999)
        self.assertGreater(nominal_likelihood((-1,), F(1), probe, 1), .999)

    def test_whole_family_and_sequential_records(self):
        belief = prior_empty()
        truth = (1, -1, 1)
        for level, sign in enumerate(truth, 1):
            probe = choose_probe(belief)
            belief = update_belief(extend_belief(belief), level, F(3, 2), probe, int(sign == -1))
            self.assertEqual(len(belief), 2 ** level)
            self.assertEqual(current_target_decision(belief), sign)
        self.assertEqual(belief.workcount, 2 + 4 + 8)
        self.assertGreater(belief[truth], .98)

    def test_nonunit_probe_nonfinite_belief_and_family_cap(self):
        p = prior_empty()
        for _ in range(3): p = extend_belief(p)
        with self.assertRaises(ValueError): extend_belief(p)
        with self.assertRaises(ValueError): extend_belief({(): float('nan')})
        with self.assertRaises(ValueError): update_belief({(1,): 1.}, 0, F(1), null_probe(()), 0)
        with self.assertRaises(ValueError): update_belief({(1,): 1.}, 1, F(1), (F(1),)*4, 0)
        with self.assertRaises(ValueError): nominal_likelihood((1,), F(1), null_probe(()), True)


if __name__ == '__main__': unittest.main()
