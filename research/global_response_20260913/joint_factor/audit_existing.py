"""Audit the existing cubic factor against the global-response H6 fixture.

This is deliberately an audit, not a new optimizer: it proves that the
published cubic artifact uses the same rational Hamiltonian and that its
lower replay constructs only exterior-power blocks.  Discovery costs remain
charged separately.
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from research.certificate_scaling.cubic_interval_replay import run

FIXTURE = ROOT / "results/molecular_collective_20260913/campaign/h6/fixture.json"
CERT = ROOT / "results/lambda_runs/cubic_precision/downloaded/precision_scs/h6_row-column/certificate.json"
REF = ROOT / "results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json"
PROOF = ROOT / "results/certificate_scaling/cubic_precision/wedge_spectral/h6_conditioned/witness.json"

def main(out: Path) -> dict:
    fixture = json.loads(FIXTURE.read_text())
    cert = json.loads(CERT.read_text())
    hhash = lambda d: hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()
    if (fixture["modes"], fixture["particles"]) != (cert["modes"], cert["particles"]):
        raise ValueError("sector mismatch")
    if hhash(fixture["hamiltonian"]) != hhash(cert["hamiltonian"]):
        raise ValueError("Hamiltonian mismatch")
    rec = run(CERT, REF, out.with_suffix(".interval.json"), method="spectral", proof=PROOF)
    witness = json.loads(PROOF.read_text())
    dims = [sum(len(c["indices"]) for c in b["components"]) for b in witness["bodies"]]
    result = {
        "kind": "existing_cubic_factor_audit_v1",
        "fixture": str(FIXTURE), "certificate": str(CERT),
        "hamiltonian_sha256": hhash(fixture["hamiltonian"]),
        "same_frozen_hamiltonian": True,
        "interval": rec["width_float"], "lower": rec["lower"], "upper": rec["upper"],
        "wedge_body_dimensions": dims,
        "many_body_states_enumerated_by_lower_replay": 0,
        "discovery_is_input_only": True,
        "discovery_enumerates_many_body_states": False,
        "discovery_scope": "full cubic SOS dictionary and numerical solve from Hamiltonian/sector inputs; expensive dictionary construction remains charged",
        "conclusion": "A complete compact-factor lower replay and determinant-free Hamiltonian-input discovery already exist; the unresolved issue is cheap small-factor discovery.",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    return result

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(); p.add_argument("--out", type=Path, required=True)
    print(json.dumps(main(p.parse_args().out), indent=2))
