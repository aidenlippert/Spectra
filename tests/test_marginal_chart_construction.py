import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_coherent_tree import CoherentCharge, build, replay
from tests.test_marginal_coherent_tree import fixture


class ChartConstructionTests(unittest.TestCase):
    def test_chart_construction_checks_leaves_and_complete_coverage(self):
        data = fixture(); data['joint_line_bound'] = True
        oracle = CoherentCharge(data)
        original = oracle.branch_lower
        visited = []
        def checked(mask, bits):
            free = [i for i in range(oracle.modes) if not mask & (1 << i)]
            self.assertTrue(not free or (len(free) == 2 and free[0] % 2 == free[1] % 2))
            visited.append((mask, bits))
            return original(mask, bits)
        with patch.object(oracle, 'branch_lower', side_effect=checked):
            tree, cost = oracle.cover(F(2), chart_construction=True)
        self.assertTrue(visited)
        self.assertEqual(cost['covered_Q_configurations'], 30)
        self.assertTrue(cost['construction_chart_only'])
        certificate = dict(data, kind='valence_coherent_tree_v1', tree=tree, max_nodes=10000, target_lower='2')
        self.assertEqual(replay(certificate)['covered_Q_configurations'], 30)
        bad = copy.deepcopy(certificate); bad['target_lower'] = '1000'
        with self.assertRaises(ValueError): replay(bad)
        with self.assertRaises(ValueError): CoherentCharge(data).cover(F(1000), chart_construction=True)
        with self.assertRaises(ValueError): CoherentCharge(data).cover(F(2), chart_construction=True, max_nodes=1)
        with self.assertRaises(ValueError): CoherentCharge(data).cover(F(2), tree, chart_construction=True)
        with self.assertRaises(ValueError): CoherentCharge(fixture()).cover(F(2), chart_construction=True)

    def test_simplex_tree_roundtrip_and_strict_selectors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'source.json'; source.write_text(json.dumps(fixture()))
            result = build(source, root/'proof', F(2), joint_simplex=True, chart_construction=True)
            c = json.loads((root/'proof/certificate.json').read_text())
            self.assertEqual(replay(c), result)
            self.assertEqual(result['covered_Q_configurations'], 30)
            self.assertGreater(result['joint_simplex_branches'], 0)
            self.assertGreaterEqual(result['joint_simplex_maximum_assignments'], 3)
            construction = json.loads((root/'proof/construction_receipt.json').read_text())
            self.assertEqual(construction['fully_conditioned_transition_sources'], result['fully_conditioned_transition_sources'])
            bad = copy.deepcopy(c); bad['joint_simplex_bound'] = 1
            with self.assertRaises(ValueError): replay(bad)
            bad = copy.deepcopy(c); bad['tree'] = ['bound']
            with self.assertRaises(ValueError): replay(bad)


if __name__ == '__main__': unittest.main()
