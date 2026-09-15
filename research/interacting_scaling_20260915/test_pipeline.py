from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
from research.interacting_scaling_20260915 import pipeline


class PipelineTests(unittest.TestCase):
    def test_deadline_prevents_new_work_and_stage_failure_stops(self):
        with tempfile.TemporaryDirectory() as temporary, patch.object(pipeline, 'OUT', Path(temporary)), patch.object(pipeline.subprocess, 'run') as run:
            result = pipeline.execute('expired', [('a', 30, ['unused'])], deadline=-1)
            self.assertFalse(result['completed'])
            run.assert_not_called()
            run.return_value = SimpleNamespace(returncode=1)
            result = pipeline.execute('failure', [('a', 30, ['first']), ('b', 30, ['second'])])
            self.assertFalse(result['completed'])
            self.assertEqual(run.call_count, 1)
            self.assertEqual(result['last_stage'], 'a')

    def test_global_screen_selection_and_continuation(self):
        with tempfile.TemporaryDirectory() as temporary:
            case = Path(temporary)
            with self.assertRaises(ValueError): pipeline.steps_for(case, construct_nonsinglet=True, direct_nonsinglet=True)
            (case/'continuation_seed.json').write_text('{}')
            steps = pipeline.steps_for(case, discover_state=True, direct_nonsinglet=True)
            self.assertIn('--initial', steps[0][2])
            labels = [label for label, _, _ in steps]
            self.assertLess(labels.index('upper'), labels.index('magnetic_initialize'))
            self.assertLess(labels.index('magnetic_attach'), labels.index('prepare'))


if __name__ == '__main__': unittest.main()
