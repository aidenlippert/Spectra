import unittest
from fractions import Fraction

from experiments.v7_grouping_probe import residual_terms, _pack, verify, run


class GroupingProbeTests(unittest.TestCase):
    def test_residual_coefficients_are_exact_and_untruncated(self):
        terms = residual_terms(3, 2)
        assert terms and all(isinstance(t.coeff, Fraction) for t in terms)
        assert all(t.pauli != 'III' for t in terms)


    def test_partition_witness_is_independently_validated(self):
        terms = residual_terms(4, 3)
        groups, comparisons = _pack(terms, terms)
        verify(groups, terms)
        assert comparisons > 0


    def test_probe_preserves_separate_costs_and_no_failures(self):
        result = run()
        assert not result['failures']
        assert all('norm_upper' in row and 'construction_comparisons' in row and 'checker_pair_checks' in row for row in result['rows'])
