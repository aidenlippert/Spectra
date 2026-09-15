import json
import unittest
from pathlib import Path
from research.composable_response_20260913.joint import load_case, check_gap

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/global_response_20260913/exceptional"

class ExceptionalTests(unittest.TestCase):
    def test_h6_and_h8_replay_have_no_many_body_states(self):
        for name in ("h6", "h8"):
            data, tail, _, _ = load_case(name)
            cert = json.loads((OUT / f"{name}_high3.json").read_text())["certificate"]
            receipt = check_gap(data, tail, cert)
            self.assertEqual(receipt["many_body_states_enumerated"], 0)
            self.assertEqual(receipt["many_body_matrix_entries"], 0)

    def test_wrong_sector_partition_is_rejected(self):
        data, tail, _, _ = load_case("h6")
        cert = json.loads((OUT / "h6_high3.json").read_text())["certificate"]
        cert["last_spatial_orbitals"] = 1
        with self.assertRaises(ValueError):
            check_gap(data, tail, cert)

if __name__ == "__main__":
    unittest.main()
