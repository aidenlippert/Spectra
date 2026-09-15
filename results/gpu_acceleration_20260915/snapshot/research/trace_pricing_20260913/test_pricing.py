"""Check the physical trace metric, proposal rounding and covariance."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

import numpy as np

from research.certificate_scaling.commutator_dual_witness import seed
from research.joint_patterns_20260913.dual import integer_grams
from research.spin_completion_20260913.discovery import DIRECTION_DENOMINATOR
from research.trace_pricing_20260913.discovery import Model, normalized_modes, pricing_basis


class TracePricing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        prior = root/'results/molecular_collective_20260913/campaign/h6'
        cls.model = Model(json.loads((prior/'fixture.json').read_text()), json.loads((prior/'rank_10/tail.json').read_text()))
        source = root/'results/spin_completion_20260913/campaign/adaptive/round_4_batch_2'
        cls.dual = json.loads((source/'dual_proposal.json').read_text()); cls.span = json.loads((source/'span.json').read_text())

    def test_cached_metric_matches_exact_fixed_N_trace(self):
        m = self.model; trace = {w: seed(w, 12, 6) for w in m.rows}
        for gid in (0, 2):
            matrices, den = integer_grams(m.frames[gid]['polynomials'], [trace], 'anticommutator')
            exact = np.array([[float(F(v, den)) for v in row] for row in matrices[0]])
            self.assertTrue(np.allclose(exact, m.trace_metrics[gid], atol=2e-14, rtol=2e-14))

    def test_generalized_spectrum_is_invariant_to_frame_rescaling(self):
        m = self.model; gid = 2; C = m.frame_coefficients[gid]; y = np.array(self.dual['values'])
        monomial_gram = (m.price_maps[gid]@y).reshape(C.shape[0], C.shape[0])
        trace = m.monomial_traces[gid]
        Q, _, T = pricing_basis(C, trace)
        eigenvalues, _ = normalized_modes(Q.T@monomial_gram@Q, T)
        S = np.diag(np.geomspace(.1, 10., C.shape[1]))
        Q2, _, T2 = pricing_basis(C@S, trace)
        changed, _ = normalized_modes(Q2.T@monomial_gram@Q2, T2)
        self.assertTrue(np.allclose(eigenvalues, changed, atol=2e-9, rtol=2e-9))

    def test_rounded_proposals_are_negative_and_new_in_operator_space(self):
        m = self.model; proposals, _ = m.price(self.dual, self.span)
        self.assertEqual([len(g) for g in proposals], [2]*8)
        for gid, group in enumerate(proposals):
            C = m.frame_coefficients[gid]
            old = [C@(np.array(e['vector'])/DIRECTION_DENOMINATOR) for e in self.span if e['group'] == gid]
            for entry in group:
                self.assertLess(entry['rounded_ordinary_moment'], 0.)
                self.assertGreater(entry['rounded_trace_size'], 0.)
                self.assertAlmostEqual(entry['proposed_moment'], entry['rounded_normalized_moment'], places=7)
                column = C@(np.array(entry['vector'])/DIRECTION_DENOMINATOR)
                before = np.linalg.matrix_rank(np.column_stack(old), tol=1e-7); old.append(column)
                self.assertGreater(np.linalg.matrix_rank(np.column_stack(old), tol=1e-7), before)

    def test_invalid_duals_and_singular_metrics_are_refused(self):
        raw = deepcopy(self.dual); raw['values'][0] = float('inf')
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            self.model.price(raw, self.span)
        raw = deepcopy(self.dual); raw['rows'] = raw['rows'][1:]
        with self.assertRaisesRegex(ValueError, 'different coefficient map'):
            self.model.price(raw, self.span)
        with self.assertRaises(np.linalg.LinAlgError):
            normalized_modes(np.eye(2), np.diag([1., 0.]))
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            normalized_modes(np.eye(2), np.diag([1., float('nan')]))
        with self.assertRaisesRegex(ValueError, 'Rank-deficient'):
            pricing_basis(np.array([[1., 1.], [0., 0.]]), np.eye(2))


if __name__ == '__main__':
    unittest.main()
