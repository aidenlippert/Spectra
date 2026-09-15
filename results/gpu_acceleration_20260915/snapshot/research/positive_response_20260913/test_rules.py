"""Focused acceptance, refusal, and independent-action tests for the new rules."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest

from research.positive_response_20260913 import block_response as block
from research.positive_response_20260913 import coercivity as gap
from research.positive_response_20260913 import molecular_diagnostic as molecular


class PositiveResponseTests(unittest.TestCase):
    def setUp(self):
        self.model = next(block.cases())
        self.cert = block.reference.propose(self.model['reference'])

    def test_nonscalar_large_error_control(self):
        receipt = block.check(self.model, self.cert)
        control = block.exact_control(self.model, self.cert)
        self.assertLess(F(receipt['width']), F(1, 10**12))
        self.assertLess(F(control['scalar_error_test_counterexample']['K_minus_EdaggerE_over_d']), 0)
        # This case has constant diagonal but a nonzero off-diagonal: it is nonscalar.
        self.assertEqual(F(control['eliminated_block_diagonal_spread']), 0)
        self.assertNotEqual(F(control['eliminated_block_offdiagonal_example']['coefficient']), 0)

    def test_unproved_bulk_is_refused(self):
        for changes in ({'lambda_eliminated': '-1/100'}, {'lambda_retained': 0.1}, {'arbitrary_bath': 'allowed'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                block.check(self.model | changes, self.cert)

    def test_gap_and_schur_failures_are_refused(self):
        for lower in ('3/10', '0'):
            with self.subTest(lower=lower), self.assertRaises(ValueError):
                block.check(self.model, self.cert | {'lower': lower})

    def test_old_checker_does_not_accept_the_new_rule(self):
        with self.assertRaises(ValueError):
            block.reference.check(self.model, self.cert)
        nested = deepcopy(self.model); nested['reference']['extra_term'] = '1'
        with self.assertRaises(ValueError):
            block.check(nested, self.cert)

    def test_zero_bulk_reduces_to_the_reference(self):
        model = self.model | {'lambda_retained': '0', 'lambda_eliminated': '0'}
        self.assertEqual(block.check(model, self.cert)['lower'], block.reference.check(model['reference'], self.cert)['lower'])


class GapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        src = gap.ROOT/'results/molecular_collective_20260913/campaign/h6'
        cls.data = json.loads((src/'fixture.json').read_text())
        cls.tail = json.loads((src/'rank_10/tail.json').read_text())
        cls.cert = json.loads((gap.OUT/'gap_certificate.json').read_text())

    def test_accepted_exact_reconstruction(self):
        receipt = gap.check(self.data, self.tail, self.cert)
        self.assertEqual(receipt['basis_dimension'], 35)
        self.assertGreater(F(receipt['double_occupancy_bulk_gap_Ha']), 0)

    def test_false_square_root_bound_is_refused(self):
        cert = deepcopy(self.cert); cert['sqrt_upper'][0] = '0'
        with self.assertRaisesRegex(ValueError, 'norm bound'):
            gap.check(self.data, self.tail, cert)

    def test_false_reconstruction_is_refused(self):
        cert = deepcopy(self.cert); cert['reconstruction'][0][0] = str(F(cert['reconstruction'][0][0])+1)
        with self.assertRaisesRegex(ValueError, 'reconstruct'):
            gap.check(self.data, self.tail, cert)

    def test_forward_bracket_reference_is_refused(self):
        cert = deepcopy(self.cert); index = next(i for i, row in enumerate(cert['nodes']) if 'bracket' in row)
        cert['nodes'][index] = {'bracket': [0, index]}
        with self.assertRaisesRegex(ValueError, 'earlier'):
            gap.check(self.data, self.tail, cert)

    def test_wrong_binding_and_nonpositive_weights_are_refused(self):
        with self.assertRaisesRegex(ValueError, 'binding'):
            gap.check(self.data, self.tail, self.cert | {'tail_sha256': '0'*64})
        tail = deepcopy(self.tail); tail['factors'][0]['weight'] = '0'
        with self.assertRaisesRegex(ValueError, 'positive'):
            gap.check(self.data, tail, self.cert)

    def test_rounding_sqrt_up_does_not_lose_tiny_positive_increment(self):
        x = F(9, 4)+F(1, 3*2**100)
        self.assertGreater(gap.sqrt_upper(x), F(3, 2))
        self.assertGreaterEqual(gap.sqrt_upper(x)**2, x)


class MolecularTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (gap.ROOT/'results/response_consistency_20260913/separator.json').read_bytes()
        cls.parent = (gap.ROOT/'results/trace_pricing_20260913/full_dual/witness.json').read_bytes()
        cls.cert = json.loads((gap.OUT/'spin_projected_separator.json').read_text())

    def test_spin_scalar_rejection_and_physical_controls(self):
        receipt = molecular.check_projected(self.cert, self.source, self.parent)
        self.assertLess(F(receipt['spin_scalar_moment']), 0)
        self.assertEqual(receipt['generated_amplitudes_in_controls'], 64)
        self.assertGreater(F(receipt['singlet_physical_control']['exact_positive_expectation']), 0)

    def test_fitted_coefficient_is_not_accepted_as_a_symmetry_projection(self):
        cert = deepcopy(self.cert); cert['coefficients'][0] += 1
        with self.assertRaisesRegex(ValueError, 'spin projection'):
            molecular.check_projected(cert, self.source, self.parent)

    def test_wrong_projection_source_hash_is_refused(self):
        cert = self.cert | {'source_separator_sha256': '0'*64}
        with self.assertRaisesRegex(ValueError, 'hash-bound'):
            molecular.check_projected(cert, self.source, self.parent)

    def test_casimir_identity_also_holds_at_odd_spatial_size(self):
        self.assertTrue(molecular.casimir_identity(3)['all_particle_sectors_CAR_identity'])

    def test_one_body_spin_asymmetry_is_refused(self):
        p = {'spatial': 2, 'one': {((1, 0), (0, 0)): F(1)}}
        with self.assertRaisesRegex(ValueError, 'spin-independent'):
            molecular.spatial_one_body(p)

    def test_molecular_gap_survives_removing_the_weak_commutator_contribution(self):
        src = gap.ROOT/'results/molecular_collective_20260913/campaign/h6'
        data = json.loads((src/'fixture.json').read_text()); tail = json.loads((src/'rank_10/tail.json').read_text())
        p = molecular.extract(data, tail['center_number']); constant, matrix = molecular.spatial_one_body(p)
        prior_upper = json.loads((gap.ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json').read_text())
        ell = F(molecular.tail_replay(data, tail)['lower_operator_shift_Ha'])
        rows = molecular.complement_bounds(constant, matrix, F(0), ell, F(prior_upper['upper'])-F(1, 625))
        self.assertEqual([r['doubly_occupied_spatial_orbital'] for r in rows if r['strict_gap_certified']], [5])
        self.assertGreater(F(rows[5]['shifted_D_lower_Ha']), F(19, 100))
        self.assertEqual(rows[5]['eliminated_sector_dimension'], 210)


class StandaloneSectorGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        src = gap.ROOT/'results/molecular_collective_20260913/campaign/h6'
        cls.data = json.loads((src/'fixture.json').read_text())
        cls.tail = json.loads((src/'rank_10/tail.json').read_text())
        cls.cert = json.loads((gap.OUT/'sector_gap_certificate.json').read_text())

    def test_independent_molecular_gap_is_accepted(self):
        receipt = molecular.sector_gap_check(self.data, self.tail, self.cert)
        self.assertGreater(F(receipt['accepted_delta_Ha']), F(19, 100))
        self.assertEqual(receipt['lie_bracket_nodes_required'], 0)

    def test_oversized_or_nonpositive_gap_is_refused(self):
        for delta in ('1', '0', '-1'):
            with self.subTest(delta=delta), self.assertRaises(ValueError):
                molecular.sector_gap_check(self.data, self.tail, self.cert | {'delta_Ha': delta})

    def test_wrong_sector_and_binding_are_refused(self):
        for changes in ({'doubly_occupied_spatial_orbital': 0}, {'doubly_occupied_spatial_orbital': 6},
                        {'fixture_sha256': '0'*64}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                molecular.sector_gap_check(self.data, self.tail, self.cert | changes)


if __name__ == '__main__':
    unittest.main()
