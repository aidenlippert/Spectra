from copy import deepcopy
from fractions import Fraction as F
from unittest.mock import patch
import unittest
from experiments.v4_planner import (new_work, h1_values, action_values_h2, eval_two_step_policy,
                                    verify_fixed_policy_certificate, verify_policy_certificate)


class PlannerTests(unittest.TestCase):
    def test_known_two_model_risks_and_work(self):
        prior, signs = (F(1, 2), F(1, 2)), (-1, 1)
        table = ((F(1, 2), F(0)), (F(1, 2), F(1)))
        work = new_work()
        self.assertEqual(h1_values(prior, table, signs, work), (F(1, 2), F(0)))
        self.assertEqual(work['likelihood_weight_products'], 8)
        self.assertEqual(action_values_h2(prior, table, signs), (F(0), F(0)))
        cert = eval_two_step_policy(prior, table, signs, 0)
        self.assertEqual(cert['risk'], '0')
        with patch('experiments.v4_planner.eval_two_step_policy', side_effect=AssertionError('generator called')):
            self.assertTrue(verify_fixed_policy_certificate(cert, prior, table, signs))
            self.assertTrue(verify_policy_certificate(cert, prior, table, signs))

    def test_impossible_branch_and_corruption(self):
        prior, signs = (F(1, 2), F(1, 2)), (-1, 1)
        table = ((F(1), F(1, 2)), (F(1), F(1, 2)))
        cert = eval_two_step_policy(prior, table, signs, 0)
        self.assertIsNone(cert['branches'][0]['second_action'])
        self.assertEqual(action_values_h2(prior, table, signs), (F(1, 2), F(1, 2)))
        self.assertTrue(verify_fixed_policy_certificate(cert, prior, table, signs))
        for field, value in [('risk', '0'), ('first_action', .5)]:
            bad = deepcopy(cert); bad[field] = value
            self.assertFalse(verify_fixed_policy_certificate(bad, prior, table, signs))
        bad = deepcopy(cert); bad['branches'][1]['probability'] = '1/2'
        self.assertFalse(verify_fixed_policy_certificate(bad, prior, table, signs))

    def test_leaf_decision_changes_actual_risk(self):
        prior, signs = (F(1, 2), F(1, 2)), (-1, 1)
        table = ((F(0),), (F(1),))
        cert = eval_two_step_policy(prior, table, signs, 0)
        cert['branches'][1]['leaf_predictions'][1] = -1
        self.assertFalse(verify_fixed_policy_certificate(cert, prior, table, signs))

    def test_strict_inputs(self):
        with self.assertRaises(ValueError): h1_values((1.0,), ((F(1, 2),),), (1,))
        with self.assertRaises(ValueError): h1_values((F(1),), ((F(2),),), (1,))
        with self.assertRaises(ValueError): eval_two_step_policy((F(1),), ((F(1),),), (1,), True)


if __name__ == '__main__': unittest.main()
