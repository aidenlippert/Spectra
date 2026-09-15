import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from research.interacting_scaling_20260915 import capacity


class CapacityTests(unittest.TestCase):
    def test_prepared_maps_and_changed_adapter_refuse_an_envelope_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); case = root/'case'; case.mkdir()
            design = {'max_Gram_entries': 2000000, 'max_coefficient_rows': 180000}
            path = case/'design.json'; path.write_text(json.dumps(design))
            protocol = root/'protocol.json'
            record = {'adapter_sha256': hashlib.sha256(Path(capacity.__file__).read_bytes()).hexdigest(),
                'new_Gram_entry_cap': 4000000, 'new_coefficient_row_cap': 400000}
            protocol.write_text(json.dumps(record))
            with patch.object(capacity, 'PROTOCOL', protocol):
                (case/'prepared').mkdir()
                with self.assertRaises(ValueError): capacity.apply(case)
                self.assertEqual(json.loads(path.read_text()), design)
                (case/'prepared').rmdir()
                record['adapter_sha256'] = '0'*64
                protocol.write_text(json.dumps(record))
                with self.assertRaises(ValueError): capacity.apply(case)
                self.assertEqual(json.loads(path.read_text()), design)


if __name__ == '__main__': unittest.main()
