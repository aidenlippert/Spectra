import copy
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_boolean_line import add, multiply, absolute, exponential, compile_line
from experiments.marginal_coherent_tree import CoherentCharge, build, replay
from tests.test_marginal_coherent_tree import fixture, direct_rows


class BooleanLineTests(unittest.TestCase):
    def test_exact_ring_absolute_and_exponential_with_sign_changes(self):
        a, b = (F(-2), F(5)), (F(3, 2), F(-7, 3))
        for x in (0, 1):
            for poly, expected in [(multiply(a, b), (a[0]+x*a[1])*(b[0]+x*b[1])),
                                   (absolute(a), abs(a[0]+x*a[1])),
                                   (exponential(F(2, 3), (F(-2), F(3))), F(2, 3)**(-2+3*x))]:
                self.assertEqual(poly[0]+x*poly[1], expected)
        for base, power in [(F(0), (F(1), F(0))), (F(-1), (F(1), F(0))),
                            (F(2), (F(1, 2), F(0)))]:
            with self.assertRaises(ValueError): exponential(base, power)

    def test_every_two_mode_chart_matches_physical_rows_in_interfering_model(self):
        o = CoherentCharge(fixture()); rows = direct_rows(o)
        checked = 0
        for u, v in combinations(range(o.modes), 2):
            if u % 2 != v % 2: continue
            mask = o.all_bits ^ ((1 << u) | (1 << v))
            for bits in {s & mask for s in rows}:
                if o.close(mask, bits) != (mask, bits): continue
                compiled = compile_line(o, mask, bits)
                self.assertIsNotNone(compiled); poly, count = compiled
                states = [bits | (1 << v), bits | (1 << u)]
                physical = [rows[s] for s in states if s in rows]
                self.assertEqual(poly[0]+min(F(0), poly[1]), min(physical))
                for x, s in enumerate(states):
                    if s in rows: self.assertEqual(poly[0]+x*poly[1], rows[s])
                self.assertEqual(count['feasible_assignments'], 2)
                checked += 1
        self.assertGreater(checked, 20)
        self.assertIsNone(compile_line(o, 0, 0))

    def test_saved_h6_frontier_is_repaired_without_determinant_actions(self):
        root = Path(__file__).resolve().parents[1]/'results/marginal_h6/coherent_charge_range2_polynomial'
        c = json.loads((root/'certificate.json').read_text()); failure = json.loads((root/'frontier_branch_obstruction.json').read_text())
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                   side_effect=AssertionError('No determinant action')):
            o = CoherentCharge(c); row, counts = compile_line(o, failure['mask'], failure['bits'])
        self.assertEqual([row[0], sum(row)], list(map(F, failure['exact_coherent_row_lowers'])))
        self.assertGreater(row[0]+min(F(0), row[1]), F(c['target_lower']))
        self.assertGreater(counts['metric_scalar_endpoint_evaluations'], 0)

    def test_complete_tree_exports_joint_operation_costs_and_refuses_false_bound(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'source.json'; source.write_text(json.dumps(fixture()))
            with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                       side_effect=AssertionError('No determinant action')):
                r = build(source, root/'proof', F(2), source_polynomial=True, joint_line=True)
                c = json.loads((root/'proof/certificate.json').read_text())
                self.assertEqual(replay(c), r)
                construction = json.loads((root/'proof/construction_receipt.json').read_text())
                self.assertEqual(construction['covered_Q_configurations'], 30)
                self.assertGreaterEqual(construction['fully_conditioned_transition_sources'], r['fully_conditioned_transition_sources'])
            self.assertGreater(r['joint_line_branches'], 0)
            self.assertEqual(r['joint_line_assignment_capacity'], 2*r['joint_line_branches'])
            self.assertEqual(r['covered_Q_configurations'], 30)
            bad = copy.deepcopy(c); bad['target_lower'] = '1000'
            with self.assertRaises(ValueError): replay(bad)
            bad = copy.deepcopy(c); bad['joint_line_bound'] = 1
            with self.assertRaises(ValueError): replay(bad)


if __name__ == '__main__': unittest.main()
