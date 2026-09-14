from fractions import Fraction as F
from itertools import combinations
import unittest
from unittest.mock import patch

from experiments.marginal_boolean_line import compile_line
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_occupation_simplex import chart_modes, compile_simplex
from tests.test_marginal_coherent_tree import fixture, direct_rows


class OccupationSimplexTests(unittest.TestCase):
    def test_all_supported_charts_against_independent_car_rows(self):
        oracle = CoherentCharge(fixture()); rows = direct_rows(oracle)
        sizes, holes, checked_p = set(), set(), 0
        for k in (2, 3, 4):
            for free in combinations(range(oracle.modes), k):
                if len({i % 2 for i in free}) != 1: continue
                mask = oracle.all_bits ^ sum(1 << i for i in free)
                for bits in {state & mask for state in rows}:
                    if oracle.close(mask, bits) != (mask, bits): continue
                    chart = chart_modes(oracle, mask, bits)
                    if chart is None: continue
                    with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                               side_effect=AssertionError('No determinant actions')):
                        values, counts = compile_simplex(oracle, mask, bits)
                    _, one_hole = chart
                    states = [bits | sum(1 << i for i in free if (i != selected if one_hole else i == selected))
                              for selected in free]
                    self.assertEqual(min(values), min(rows[s] for s in states if s in rows))
                    for value, state in zip(values, states):
                        if state in rows: self.assertEqual(value, rows[state])
                        else:
                            checked_p += 1
                            self.assertEqual(value, sum(map(abs, oracle.oracle.h.values()), F(0)))
                    self.assertEqual(counts['feasible_assignments'], k)
                    self.assertEqual(counts['amplitude_scalar_endpoint_evaluations'], k*counts['transition_groups_used'])
                    self.assertEqual(counts['metric_scalar_endpoint_evaluations'], k*counts['metric_factors_touched'])
                    if k == 2:
                        line, _ = compile_line(oracle, mask, bits)
                        self.assertEqual(values, (sum(line), line[0]))
                    sizes.add(k); holes.add(one_hole)
        self.assertEqual(sizes, {2, 3}); self.assertEqual(holes, {False, True})
        self.assertGreater(checked_p, 0)

    def test_ineligible_and_invalid_chart_refusal(self):
        o = CoherentCharge(fixture())
        self.assertIsNone(compile_simplex(o, 0, 0))
        # Two free modes of opposite spin are not a simplex chart.
        mask = o.all_bits ^ 3
        self.assertIsNone(compile_simplex(o, mask, (1 << 2) | (1 << 3)))
        # Same-spin one-electron chart, but the other spin is underfilled.
        mask = o.all_bits ^ ((1 << 0) | (1 << 2))
        with self.assertRaises(ValueError): compile_simplex(o, mask, 1 << 4)


if __name__ == '__main__': unittest.main()
