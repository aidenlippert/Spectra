"""Boundary and cross-implementation checks for the exploratory policy run."""
from fractions import Fraction as F
from copy import deepcopy
import unittest
from unittest.mock import patch
import numpy as np
from experiments.v4_run import design_case, features, context
from experiments.v4_moments import signed_gram_plan
from experiments.v4_planner import (compile_fixed_actions, action_values_h2,
                                    verify_fixed_policy_certificate, new_work)


class RunTests(unittest.TestCase):
    def test_hidden_selected_model_does_not_enter_features(self):
        task = context(7001, 3, True, 101)
        before, _ = features(task)
        changed = deepcopy(task)
        changed['evaluator_hidden_model'] = 7 - task['evaluator_hidden_model']
        changed['seed'] = -1
        changed['history'] = [{'arbitrary': 'not used by the predictor'}]
        after, _ = features(changed)
        np.testing.assert_array_equal(before, after)

    def test_features_and_source_do_not_call_teacher(self):
        with patch('experiments.v4_run.signed_gram_plan', side_effect=AssertionError('teacher access')):
            task = context(7002, 4, False, 103)
            values, work = features(task)
        self.assertEqual(values.shape, (len(task['actions']), 12))
        self.assertGreater(work['likelihood_weight_products'], 0)

    def test_direct_gram_certificate_no_child_search(self):
        for seed in range(7101, 7105):
            case = design_case(seed, 2 + seed % 3, bool(seed % 2), 103)
            task, plan = case['task'], case['teacher']
            p, L, signs = task['prior'], task['table'], task['signs']
            self.assertEqual(plan['values'], action_values_h2(p, L, signs))
            work = new_work()
            a = plan['first_action']
            with patch('experiments.v4_planner._h1', side_effect=AssertionError('redundant search')):
                cert = compile_fixed_actions(p, L, signs, a, plan['child_actions'][a], work)
                self.assertTrue(verify_fixed_policy_certificate(cert, p, L, signs))
            self.assertEqual(F(cert['risk']), min(plan['values']))
            self.assertEqual(work['one_step_action_scores'], 0)
            self.assertEqual(work['posterior_normalizations'], 0)
            self.assertLessEqual(work['likelihood_weight_products'], 6 * len(p))
            bad = deepcopy(cert); bad['risk'] = str(F(cert['risk']) + F(1, 100))
            self.assertFalse(verify_fixed_policy_certificate(bad, p, L, signs))

    def test_zero_branch_certificate_and_invalid_actions(self):
        p, L, signs = (F(1), F(0)), ((F(1),), (F(0),)), (1, -1)
        cert = compile_fixed_actions(p, L, signs, 0, (0, 0))
        self.assertEqual(cert['branches'][0]['second_action'], None)
        self.assertTrue(verify_fixed_policy_certificate(cert, p, L, signs))
        with self.assertRaises(ValueError): compile_fixed_actions(p, L, signs, 0, (True, 0))


if __name__ == '__main__': unittest.main()
