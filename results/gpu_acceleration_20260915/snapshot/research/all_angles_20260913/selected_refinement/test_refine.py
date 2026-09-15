from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from research.all_angles_20260913.selected_refinement.refine import CompiledHamiltonian, ROOT, run
from research.certificate_scaling.streaming_reference_upper import upper


class RefinementTests(unittest.TestCase):
    def test_compilation_matches_independent_car_action(self):
        fixture = json.loads((ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json').read_text())
        compiled = CompiledHamiltonian(fixture)
        for state in (15, 51, 85, 170):
            expected = compiled.original.action(state)
            actual = compiled.action(state)
            self.assertEqual(set(actual), set(expected))
            for target, value in expected.items():
                self.assertEqual(actual[target], float(value))

    def test_rational_witness_and_refusal(self):
        fixture_path = ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json'
        fixture = json.loads(fixture_path.read_text())
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'run'
            result = run(fixture_path, out, [8,16], seconds=20)
            witness = json.loads((out/result['rows'][-1]['witness_file']).read_text())['independent_upper']
            value, _ = upper(fixture, witness)
            self.assertEqual(value, Fraction(result['rows'][-1]['upper']))
            witness['states'][0] = 0
            with self.assertRaises(ValueError):
                upper(fixture, witness)
        with self.assertRaises(ValueError):
            run(fixture_path, '/unused', [8,4])


if __name__ == '__main__':
    unittest.main()
