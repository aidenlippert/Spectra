"""Verify exact separation and reject a bound-but-nonseparating direction."""
from copy import deepcopy
from pathlib import Path
import json
import tempfile
import unittest

from research.spin_completion_20260913.diagnostic import check_separator
from research.spin_completion_20260913.raw_control import run as run_full_control


class Separators(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        prior = root/'results/molecular_collective_20260913/campaign/h6'
        cls.data = json.loads((prior/'fixture.json').read_text()); cls.tail = json.loads((prior/'rank_10/tail.json').read_text())
        cls.dual = json.loads((root/'results/spin_subspace_20260913/full_dual/witness.json').read_text())
        cls.cert = json.loads((root/'results/spin_completion_20260913/diagnostic/separator_2_0.json').read_text())

    def test_saved_new_direction_is_exactly_negative(self):
        receipt = check_separator(self.data, self.tail, self.dual, self.cert)
        self.assertLess(receipt['normalized_expectation_float'], -0.017)
        self.assertLessEqual(receipt['max_degree'], 4)

    def test_wrong_binding_and_positive_linear_condition_are_refused(self):
        cert = deepcopy(self.cert); cert['dual_sha256'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'binding failed'):
            check_separator(self.data, self.tail, self.dual, cert)
        cert = deepcopy(self.cert); cert['generators'] = [[-1, -1, -1, 0]]; cert['vector'] = [1]
        with self.assertRaisesRegex(ValueError, 'does not separate'):
            check_separator(self.data, self.tail, self.dual, cert)

    def test_full_control_preserves_existing_results_and_sidecars(self):
        for existing in ('directory', 'log', 'watchdog'):
            with self.subTest(existing=existing), tempfile.TemporaryDirectory() as directory:
                out = Path(directory)/'control'
                if existing == 'directory':
                    out.mkdir(); saved = out/'receipt.json'
                elif existing == 'log':
                    saved = out.with_suffix('.log')
                else:
                    saved = out.with_name(out.name+'_watchdog.json')
                saved.write_text('preserve this evidence')
                with self.assertRaisesRegex(FileExistsError, 'Refusing to overwrite'):
                    run_full_control(out)
                self.assertEqual(saved.read_text(), 'preserve this evidence')


if __name__ == '__main__':
    unittest.main()
