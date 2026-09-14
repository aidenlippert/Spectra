from fractions import Fraction as F
import unittest
import numpy as np
from experiments.v2_revision import Law, acquire, diagnose, sample_callback, run


class RevisionTests(unittest.TestCase):
    def test_sequential_discoveries_with_heldout_and_tied_baseline(self):
        result = run()
        for stage in result['stages']:
            self.assertTrue(stage['evaluator_exact_law_recovery'])
            self.assertTrue(stage['acquisition']['conventional_coefficients_tie'])
            self.assertEqual(stage['diagnostic']['status'], 'no_heldout_disagreement')
            fit = {tuple(r['context']) for r in stage['acquisition']['records']}
            diagnostic = {tuple(r['context']) for r in stage['diagnostic']['records']}
            self.assertFalse(fit & diagnostic)
        self.assertEqual(result['stages'][0]['new_interactions'], [(0, 1)])
        self.assertEqual(result['stages'][1]['new_interactions'], [(1, 2)])
        self.assertEqual(result['costs']['warm_shot_reduction'], 0)

    def test_cubic_failure_triggers_abstention(self):
        negative = run()['cubic_negative_control']
        self.assertTrue(negative['training_indistinguishable_from_quadratic'])
        self.assertEqual(negative['diagnostic']['status'], 'model_mismatch_abstain')
        self.assertTrue(negative['diagnostic']['disagreements'])

    def test_constant_term_and_source_only_access(self):
        law = Law((1, 0, 1), ((0, 2),), constant=1)
        rng = np.random.default_rng(41)
        calls = []
        base = sample_callback(law, F(0), rng)
        def source(x, n):
            calls.append(x)
            return base(x, n)
        acq = acquire(source, 3, eta=F(0))
        self.assertEqual(acq['model'], law)
        self.assertTrue(all(sum(x) <= 2 for x in calls))
        self.assertEqual(diagnose(source, acq, 3, eta=F(0))['status'], 'no_heldout_disagreement')

    def test_insufficient_samples_refuse_before_source_access(self):
        def forbidden(*args):
            self.fail('source accessed despite insufficient authorized shot budget')
        acq = acquire(forbidden, 4, shots=1)
        self.assertEqual(acq['status'], 'insufficient_shots')
        self.assertIsNone(acq['model'])
        with self.assertRaises(ValueError): acquire(forbidden, 9)
        with self.assertRaises(ValueError): acquire(forbidden, 4, shots=2)


if __name__ == '__main__': unittest.main()
