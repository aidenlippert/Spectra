import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from research.interacting_scaling_20260915 import cases, main_refinement


class ContinuationInputTests(unittest.TestCase):
    def test_verified_external_screen_can_initialize_a_base_that_has_none(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root/'source'; source.mkdir()
            (source/'mps').mkdir()
            (source/'fixture.json').write_text(json.dumps({'modes': 8}))
            (source/'upper.json').write_text('{}')
            (source/'mps/state.json').write_text('{}')
            # This tests input plumbing, not acceptance of a physical proof.
            external = root/'separately_checked_screen.json'
            external.write_text('{"input_plumbing_test": true}')
            with patch.object(cases, 'OUT', root/'output'):
                target = cases.initialize('new_case', source, [2], nonsinglet_source=external)
            self.assertFalse((source/'nonsinglet.json').exists())
            self.assertEqual((target/'nonsinglet.json').read_bytes(), external.read_bytes())
            receipt = json.loads((target/'input_dependencies.json').read_text())
            self.assertEqual(receipt['inputs']['nonsinglet.json']['sha256'], hashlib.sha256(external.read_bytes()).hexdigest())

    def test_main_continuation_rejects_a_changed_state_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root/'source'; source.mkdir()
            (source/'prepared').mkdir(); (source/'mps').mkdir()
            (source/'fixture.json').write_text('{}')
            (source/'mps/state.json').write_text('{"changed": true}')
            (source/'design.json').write_text('{}')
            (source/'prepared/frame.json').write_text(json.dumps({'magnetization': 0,
                'fixture_sha256': hashlib.sha256(b'{}').hexdigest(),
                'state_sha256': hashlib.sha256(b'{}').hexdigest()}))
            target = root/'target'
            with patch.object(main_refinement, 'SOURCE', source), patch.object(main_refinement, 'TARGET', target):
                with self.assertRaisesRegex(ValueError, 'different source'):
                    main_refinement.initialize()
            self.assertFalse(target.exists())


if __name__ == '__main__': unittest.main()
