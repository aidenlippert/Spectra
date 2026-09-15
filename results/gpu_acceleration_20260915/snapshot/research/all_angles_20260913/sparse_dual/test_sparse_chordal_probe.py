import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_symbolic import decode, verify


def test_original_h_export_and_malformed_refusal():
    fixture = json.loads((ROOT / "results/certificate_scaling/active_space_ladder/h4/fixture.json").read_text())
    exported = json.loads((ROOT / "results/all_angles_20260913/sparse_dual/h4_sparse/certificate.json").read_text())
    original = decode(fixture["hamiltonian"], fixture["modes"], 4)
    candidate = decode(exported["hamiltonian"], exported["modes"], 4)
    assert original == candidate
    assert verify(exported)["mode_count"] == fixture["modes"]
    malformed = copy.deepcopy(exported)
    malformed["blocks"][0]["factor"][0][0] = "corrupt"
    try:
        verify(malformed)
    except ValueError:
        pass
    else:
        raise AssertionError("malformed factor was accepted")
