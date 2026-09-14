import unittest
from fractions import Fraction
import numpy as np
from experiments.discovery import dictionary, probe_bank, design, observe, Mechanism, fit_model, parameter_intervals, benchmark


class DiscoveryTests(unittest.TestCase):
    def test_sensor_matches_known_commutator(self):
        self.assertEqual(design([('Y', 'Z'), ('Z', 'Y'), ('X', 'X')], ('X',)).tolist(), [[-2.0], [2.0], [0.0]])
        self.assertEqual(len(dictionary()), 36)

    def test_exact_parameter_bounds_and_misspecification(self):
        d = dictionary(); rows = probe_bank(); a = design(rows, d)
        truth = Mechanism(('XII', 'ZXI'), (Fraction(1, 3), Fraction(-7, 8)))
        noise = Fraction(1, 500)
        values = observe(truth, rows, np.random.default_rng(8), noise)
        intervals, status = parameter_intervals(a, values, noise)
        self.assertEqual(status, 'bounded_within_declared_family')
        for p, (lo, hi) in zip(d, intervals):
            self.assertLessEqual(lo, truth.hamiltonian().get(p, 0))
            self.assertGreaterEqual(hi, truth.hamiltonian().get(p, 0))
        outside = observe(Mechanism(('XXX',), (Fraction(1),)), rows, np.random.default_rng(2), noise)
        self.assertIsNone(parameter_intervals(a, outside, noise)[0])

    def test_unidentifiability_and_frozen_no_refit(self):
        model, status = fit_model(np.zeros((3, 2)), np.zeros(3), ('X', 'Z'), 'from_scratch')
        self.assertTrue(status['abstained'])
        model, _ = fit_model(np.eye(2), np.array([2., 3.]), ('X', 'Z'), 'frozen', {'X': 0.5})
        self.assertEqual(model, {'X': 0.5})

    def test_three_stages_two_new_mechanisms_and_transfer(self):
        result = benchmark(seed=11)
        self.assertEqual(len(result['stages']), 3)
        for i, stage in enumerate(result['stages']):
            self.assertEqual(len(stage['contexts']), 2)
            for context in stage['contexts']:
                seq = context['methods']['sequential']
                self.assertFalse(seq['abstained'])
                self.assertTrue(seq['oracle_support_correct'])
                self.assertEqual(len(seq['support']), i + 2)
                self.assertLess(seq['holdout_max_error'], .006)
                self.assertIsNotNone(context['conditional_model_norm_error'])
        self.assertTrue(result['unidentifiable_abstained'])
        self.assertEqual(result['misspecified_status'], 'out_of_family_or_noise_bound_violated')
        self.assertFalse(result['compounding_demonstrated'])


if __name__ == '__main__':
    unittest.main()
