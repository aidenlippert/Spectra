"""Spin decomposition and false-separator refusal."""
from fractions import Fraction as F
import unittest

from experiments.marginal_symbolic import add
from research.molecular_collective_20260913.core import digest
from research.joint_patterns_20260913.core import generator
from research.joint_patterns_20260913.spin_diagnostic import check_separator, spin_generator
from research.joint_patterns_20260913 import test_core as fixtures


class SpinComponents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixtures.JointPatterns.setUpClass(); cls.fixture = fixtures.JointPatterns

    def test_two_spin_components_reconstruct_each_learned_generator(self):
        patterns = self.fixture.patterns
        for k in range(len(patterns)):
            for mode in range(8):
                self.assertEqual(add(spin_generator(patterns, [k, 0, mode], 8),
                    spin_generator(patterns, [k, 1, mode], 8)), generator(patterns, [k, mode], 8))

    def test_invalid_spin_and_nonnegative_direction_refused(self):
        case = self.fixture
        with self.assertRaises(ValueError):
            spin_generator(case.patterns, [0, 2, 0], 8)
        dual = {'moments': [{'word': [], 'value': '1'}]}
        cert = {'kind': 'spin_pattern_separator_v1', 'fixture_sha256': digest(case.data),
            'tail_sha256': digest(case.tail), 'dual_sha256': digest(dual),
            'generators': [[-1, -1, 0]], 'vector': [1]}
        # {a_0^dagger,a_0}=I has expectation one and cannot be a separator.
        with self.assertRaises(ValueError):
            check_separator(case.data, case.tail, dual, cert)


if __name__ == '__main__':
    unittest.main()
