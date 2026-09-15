import unittest
from fractions import Fraction as F
from copy import deepcopy
from experiments.v5_core import Contract, LearnedPrefix, acquire_next, make_source, replay_records
from experiments.v5_complementarity import certificate, retained_policy, run as retention_run
from experiments.v5_run import run_world, world, METHODS
from experiments.v5_frame import C, S, Q, frame, dot


class ResearchTests(unittest.TestCase):
    def test_B_only_gram_is_derived_from_actual_mode_geometry(self):
        # Independently derive the residual matrix for both known B values.
        for b in (-1, 1):
            coords=[tuple(dot(frame((a,b))[1],q) for q in Q[:3]) for a in (-1,1)]
            gram=tuple(tuple(sum(v[i]*v[j] for v in coords) for j in range(3)) for i in range(3))
            self.assertEqual(gram[0],(2*S**4,F(0),F(0)))
            self.assertEqual(gram[1][1],2*S*S*C*C)
            self.assertEqual(gram[2][2],2*C*C)
            self.assertEqual(gram[1][2],-2*b*S*C*C)
            self.assertEqual(gram[1][1]*gram[2][2]-gram[1][2]**2,0)
            self.assertEqual(gram[1][1]+gram[2][2],2*(1-S**4))

    def test_exact_preparation_and_gaussian_tail_bounds(self):
        q = F(1, 2**60)
        kl_upper = q/2 + q*q/(1-q)
        self.assertLess(kl_upper, q)
        self.assertLess(F(1, 2**30), F(1, 10**9))
        null_variance = 1 + F(1, 10**18)*2*(10**14+10**10+10**6)
        self.assertLess(null_variance, F(1001,1000))
        self.assertLess(F(1001,1000)/F(19999,1000)**2, F(1,399))
        self.assertLess(F(20001,1000)/(F(19,20)*1000), F(11,500))
        self.assertEqual(C*C-F(9,16)*S*S, 0)
        self.assertGreater(S*S-F(9,16)*C*C, 0)

    def test_two_transfers_and_fair_stateful_ties(self):
        data = world(41001)
        result = {method: run_world(data, method) for method in METHODS}
        expected = {'cumulative':21,'frozen_after_first':37,'scratch':56,
                    'exact_lookup':56,'structured_retrieval':21,'stateful_bayes':21,'fixed_one_read':20}
        for method, count in expected.items():
            self.assertEqual(result[method]['total_scalar_reads'], count)
            self.assertEqual(sum(len(r['records']) for r in result[method]['source_records']), count)
        self.assertEqual([r['estimate'] for r in result['cumulative']['tasks']],
                         [r['estimate'] for r in result['structured_retrieval']['tasks']])
        self.assertEqual(result['exact_lookup']['exact_cache_hits'], 0)

    def test_corrupted_evidence_and_undeclared_actuator_mismatch(self):
        first = acquire_next(lambda *args: 0, Contract(('a',), F(1)))
        record = deepcopy(first.records[0]); record['probe'] = tuple(reversed(record['probe']))
        # A different exact unit vector must not authenticate as the performed probe.
        with self.assertRaises(ValueError): replay_records((record,))
        nominal_errors = bad_errors = 0
        for seed in range(301,333):
            contract = Contract(('a',), F(1))
            nominal, _, _ = make_source((1,), contract, seed)
            bad, ledger, _ = make_source((1,), contract, seed, actual_actuator=F(1,1000),
                                         allow_contract_violation=True)
            nominal_errors += int(acquire_next(nominal,contract).signs != (1,))
            bad_errors += int(acquire_next(bad,contract).signs != (1,))
            self.assertTrue(ledger['intentional_contract_violation'])
        self.assertEqual(nominal_errors, 0)
        self.assertGreaterEqual(bad_errors,24)

    def test_common_baseline_complementarity_and_retention_boundary(self):
        cert=certificate()
        self.assertGreater(F(cert['acquired_complementarity_lower']), F(77,1000))
        # B-only access must not silently carry A's sign or its full mode vector.
        with self.assertRaises(ValueError):
            retained_policy(lambda *args:self.fail('invalid request queried source'),
                            Contract(('a','b','c'),F(1)), 'B', (1,-1))
        study=retention_run(worlds=2,targets=2)
        self.assertEqual(study['acquisition_reads'],4)
        self.assertEqual(study['counterfactual_target_reads'],16)
        for row in study['rows']:
            for target in row['targets']:
                for ledger in target['ledgers'].values():
                    self.assertEqual(ledger['scalar_readouts'],1)


if __name__ == '__main__': unittest.main()
