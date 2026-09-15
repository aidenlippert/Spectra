"""Tests of the inferential and refusal boundaries, not only happy-path output."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import product
import unittest
import numpy as np
from experiments.v2_run import acquire, compiled_b_error, contexts, conventional_fit, run
from experiments.v2_laws import CalibrationRecord, learn, predict, verify_certificate
from experiments.v2_physics import World


class TransferTests(unittest.TestCase):
    def test_compiled_channel_matches_born_distribution(self):
        for v in (F(1, 7), F(1, 2), F(1)):
            eta, coin = compiled_b_error(v)
            for bit in (0, 1):
                plus = (1 + v * (bit == 0)) / 2
                label_one = plus * coin + 1 - plus
                self.assertEqual(label_one, eta if bit == 0 else 1 - eta)

    def test_no_truth_needed_by_acquisition_and_holdout_transfer(self):
        rng = np.random.default_rng(44)
        masks = ((1, 0, 1), (0, 1, 1))
        world = World(*masks)
        calls = []
        def a_sensor(x, n):
            calls.append(('A', x, n))
            return world.calibration('a', x, n, rng)
        def reference(x, z, action, n):
            calls.append(('reference', x, z, action, n))
            return world.target(x, z, 1, action, n, rng)
        acq = acquire(a_sensor, reference, 3, F(1, 2), F(1, 1000), F(1, 1000), rng)
        self.assertEqual((acq['a'].mask, acq['b'].mask), masks)
        self.assertEqual(acq['conventional'], masks)
        self.assertEqual([c[0] for c in calls], ['A'] * 3 + ['reference'] * 3)
        training, heldout = contexts(3)
        self.assertNotIn(acq['reference_x'], training)
        for x in heldout:
            for i, name in enumerate(('a', 'b')):
                self.assertEqual(predict(acq[name].mask, x), predict(masks[i], x))

    def test_certificate_cannot_establish_class_membership(self):
        # The nonlinear alternative agrees on every supplied calibration context.
        train, _ = contexts(2)
        records = tuple(CalibrationRecord(x, (0,)) for x in train)
        law = learn(records, F(0), F(1, 100))
        self.assertTrue(verify_certificate(records, law, F(0), F(1, 100)))
        x = (1, 1)
        nonlinear_truth = x[0] & x[1]
        self.assertNotEqual(predict(law.mask, x), nonlinear_truth)

    def test_abstain_on_insufficient_or_rank_deficient_calibration(self):
        few = (CalibrationRecord((1, 0), (0,)), CalibrationRecord((0, 1), (1,)))
        self.assertFalse(learn(few, F(2, 5), F(1, 1000)).certificate_valid)
        rank_bad = (CalibrationRecord((1, 0), (0,)), CalibrationRecord((1, 0), (0,)))
        law = learn(rank_bad, F(0), F(1, 1000))
        self.assertIsNone(law.mask)
        self.assertFalse(verify_certificate(rank_bad, law, F(0), F(1, 1000)))

    def test_wrong_prefix_breaks_B_channel(self):
        # Wrong A makes every raw reference measurement a fair coin, independent
        # of B; the channel certificate therefore cannot be unconditional.
        eta, coin = compiled_b_error(F(1, 2))
        wrong_prefix_label_one = F(1, 2) + coin / 2
        self.assertEqual(wrong_prefix_label_one, 1 - eta)
        self.assertNotEqual(wrong_prefix_label_one, eta)

    def test_small_end_to_end_result(self):
        p = {'dimension': 3, 'world_seeds': [17], 'visibility': '1/2',
             'a_channel_error': '1/10', 'delta_a': '1/1000',
             'delta_b_conditional_on_a': '1/1000', 'evaluation_trials_per_world': 128,
             'evaluation_family_failure': '1/100'}
        result = run(p)
        self.assertEqual(result['theorem']['ideal_complementarity'], '1/16')
        self.assertEqual(result['theorem']['optimal_reset_four_shot_risk'], '103/256')
        w = result['worlds'][0]
        self.assertTrue(w['both_masks_recovered'])
        self.assertTrue(w['evaluation']['conventional_action_tie'])
        self.assertEqual(w['evaluation']['exact_lookup_hits'], 0)
        self.assertEqual(w['evaluation']['exact_complementarity'], '1/16')
        self.assertLess(F(w['guarantees']['AB_unconditional_error_upper_bound']), F(3, 10))
        self.assertFalse(result['physical_experiments_performed'])


if __name__ == '__main__':
    unittest.main()
