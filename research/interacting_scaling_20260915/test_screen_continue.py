import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from research.interacting_scaling_20260915 import screen_continue as continue_screen


class MagneticContinuationTests(unittest.TestCase):
    def test_domain_and_model_mismatch_refuse_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root/'source'; source.mkdir()
            (source/'prepared').mkdir()
            (source/'fixture.json').write_text('{}')
            (source/'design.json').write_text(json.dumps({'magnetization': 1}))
            frame = {'magnetization': 0, 'fixture_sha256': hashlib.sha256(b'{}').hexdigest()}
            path = source/'prepared/frame.json'
            with patch.object(continue_screen, 'OUT', root/'output'):
                path.write_text(json.dumps(frame))
                with self.assertRaisesRegex(ValueError, 'M_S=1'):
                    continue_screen.initialize(source, 'wrong_sector')
                frame.update(magnetization=1, fixture_sha256='0'*64)
                path.write_text(json.dumps(frame))
                with self.assertRaisesRegex(ValueError, 'different Hamiltonian'):
                    continue_screen.initialize(source, 'wrong_model')
                self.assertFalse((root/'output').exists())

    def test_invalid_unbounded_run_is_refused_without_loading_a_checkpoint(self):
        with self.assertRaises(ValueError): continue_screen.run(Path('/nonexistent'), 'unbounded', 10000)


if __name__ == '__main__': unittest.main()
