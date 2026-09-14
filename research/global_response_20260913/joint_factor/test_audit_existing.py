import json
from pathlib import Path
from .audit_existing import main

def test_existing_factor_matches_fixture(tmp_path):
    r = main(tmp_path / "audit.json")
    assert r["same_frozen_hamiltonian"]
    assert r["many_body_states_enumerated_by_lower_replay"] == 0
    assert r["interval"] < 0.0016
    assert r["wedge_body_dimensions"] == [12, 66, 220]
