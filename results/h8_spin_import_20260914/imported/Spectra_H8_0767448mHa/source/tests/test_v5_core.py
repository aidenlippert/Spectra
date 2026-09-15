import unittest
from fractions import Fraction as F
from dataclasses import replace
from itertools import product
from unittest.mock import patch
import numpy as np
from experiments.v5_core import (Contract, LearnedPrefix, acquire_next, make_source,
    replay_records, perturb_probe, resource_lower_bounds, ACTUATOR_BOUND, PREPARATION_TIME)
from experiments.v5_frame import frame, null_probe


class CoreTests(unittest.TestCase):
    def test_acquisition_compiles_only_raw_bit_records(self):
        learned = LearnedPrefix()
        for level, sign in enumerate((1, -1, 1), 1):
            contract = Contract(tuple('abc'[:level]), F(3, 2))
            with patch('experiments.v5_core.frame', side_effect=AssertionError('hidden frame access')):
                learned = acquire_next(lambda n, x, value=sign: int(value == -1), contract, learned)
        expected, residual = frame((1, -1, 1))
        self.assertEqual(learned.modes, expected); self.assertEqual(learned.residual, residual)
        self.assertEqual(replay_records(learned.records), learned)
        self.assertLess(learned.failure_bound, F(1, 20))

    def test_simulated_all_eight_worlds_and_resource_ledger(self):
        for index, truth in enumerate(product((-1, 1), repeat=3)):
            contract = Contract(('a', 'b', 'c'), F(1))
            read, ledger, records = make_source(truth, contract, 31001 + index)
            learned = LearnedPrefix()
            for level in range(1, 4):
                learned = acquire_next(read, Contract(contract.lineage[:level], contract.alpha), learned)
            self.assertEqual(learned.signs, truth)
            self.assertEqual(ledger['scalar_readouts'], 3)
            self.assertEqual(ledger['source_coordinate_preparations'], 12)
            self.assertEqual(ledger['relaxation_time_units'], 3 * PREPARATION_TIME)
            self.assertEqual(len(records), 3)

    def test_actuator_bound_and_exact_positive_net_resource_bound(self):
        for length in range(3):
            for prefix in product((-1, 1), repeat=length):
                x = null_probe(prefix); applied = perturb_probe(x)
                self.assertAlmostEqual(float(np.linalg.norm(applied)), 1.)
                self.assertLessEqual(float(np.linalg.norm(applied - np.asarray(x, dtype=float))), float(ACTUATOR_BOUND))
        bounds = resource_lower_bounds()
        self.assertGreater(bounds['one_read_cold_error_lower_bound'], F(6, 25))
        self.assertGreater(bounds['adaptive_expected_cold_reads_lower_bound'], F(5, 4))
        self.assertGreater(bounds['net_expected_read_saving_per_four_target_batch_lower_bound'], 0)

    def test_refusal_before_experiment_for_precision_lineage_or_confidence(self):
        read = lambda *args: self.fail('refused request used a physical read')
        with self.assertRaises(ValueError): acquire_next(read, Contract(('a',), F(1), actuator_bound=F(1, 1000)))
        with self.assertRaises(ValueError): acquire_next(read, Contract(('a','b'), F(1)))
        with self.assertRaises(ValueError): acquire_next(read, Contract(('a',), F(1)), error_budget=F(1, 100000))
        with self.assertRaises(ValueError): acquire_next(read, Contract(('a',), F(3)))
        corrupted = replace(LearnedPrefix(), residual=(F(1),)*4)
        with self.assertRaises(ValueError): acquire_next(read, Contract(('a',), F(1)), corrupted)


if __name__ == '__main__': unittest.main()
