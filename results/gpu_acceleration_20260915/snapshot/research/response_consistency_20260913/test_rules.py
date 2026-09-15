"""Focused soundness boundaries for the two concrete proof steps."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import unittest

from research.response_consistency_20260913.response import OUT, ROOT, check as response_check, propose, cases
from research.response_consistency_20260913.separator import check as separator_check, operator, act


class ProofRules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.separator = json.loads((OUT/'separator.json').read_text())
        cls.parent = (ROOT/'results/trace_pricing_20260913/full_dual/witness.json').read_bytes()

    def test_exact_small_solvable_case(self):
        model = {'kind': 'homogeneous_central_spin_v1', 'bath_spins': 2, 'excitations': 1,
                 'epsilon': '1', 'omega': '0', 'chi': '0', 'g': '1'}
        receipt = response_check(model, {'lower': '-1', 'trial_t': '-1/2'})
        self.assertEqual(receipt['lower'], receipt['upper'])
        interacting = {**model, 'excitations': 2, 'omega': '1', 'chi': '1'}
        receipt = response_check(interacting, {'lower': '1', 'trial_t': '-1/2'})
        self.assertEqual(receipt['lower'], receipt['upper'])

    def test_response_refuses_gap_or_remainder_failure(self):
        model = cases()[0]
        for lower in ('3/10', '0'):
            with self.assertRaises(ValueError):
                response_check(model, {'lower': lower, 'trial_t': '-1'})

    def test_response_refuses_another_hamiltonian(self):
        model = cases()[0]
        with self.assertRaises(ValueError):
            response_check({**model, 'site_couplings': ['1', '2']}, propose(model))

    def test_separator_negative_and_hash_bound(self):
        result = separator_check(self.separator, self.parent)
        self.assertLess(F(result['exact_moment']), 0)
        with self.assertRaises(ValueError):
            separator_check(self.separator, self.parent+b' ')

    def test_no_invented_higher_moment(self):
        parent = json.loads(self.parent)
        parent['moments'].append({'word': [[1, 0], [1, 1], [1, 2], [0, 0], [0, 1], [0, 2]], 'value': '0'})
        blob = json.dumps(parent).encode(); cert = deepcopy(self.separator)
        cert['parent_witness_sha256'] = hashlib.sha256(blob).hexdigest()
        with self.assertRaisesRegex(ValueError, 'higher'):
            separator_check(cert, blob)

    def test_diagonal_T1_occupancy_identity(self):
        cert = {**self.separator, 'triples': [[0, 1, 2]], 'coefficients': [1], 'denominator': 1,
                'support_spatial': [0, 1]}
        _, P = operator(cert)
        for state in range(8):
            # All three empty or all three occupied: one. Otherwise: zero.
            self.assertEqual(act(P, {state: F(1)}).get(state, F(0)), int(state in (0, 7)))


if __name__ == '__main__':
    unittest.main()
