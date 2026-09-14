"""Audit compact response-lifted variational witnesses.

This is deliberately an audit of existing exact replay receipts: it does not
reuse their amplitudes as a new discovery procedure.  It records the support
and the inherited-source dependency explicitly.
"""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "compact_response_20260913"
OUT = ROOT / "results" / "global_response_20260913" / "upper_model"

def row(case, trial):
    d = trial["cases"][case]
    disc = d["trial_discovery"]
    replay = d["replay"]
    action = disc["action_cost"]
    sector = 924
    return {
        "case": d["case"],
        "sector_dimension": sector,
        "source_support": disc["initial_retained_states"],
        "response_support": disc["response_direction_states"],
        "final_support": disc["final_trial_states"],
        "source_support_fraction": disc["initial_retained_states"] / sector,
        "final_support_fraction": disc["final_trial_states"] / sector,
        "response_actions": disc["D_actions"],
        "hamiltonian_action_calls": action["Hamiltonian_action_calls"],
        "word_state_checks": action["word_state_checks"],
        "candidate_upper_Ha": replay["candidate_upper_Ha"],
        "accepted_upper_Ha": replay["accepted_upper_Ha"],
        "candidate_gain_mHa": replay["candidate_gain_mHa"],
        "source_provenance": d["source_upper_provenance"],
        "source_dependency_is_new": False,
        "scope": "Exact variational replay audit; not an enumeration-free upper-bound construction.",
    }

def main():
    trial = json.loads((SRC / "trial_discovery.json").read_text())
    rows = [row(i, trial) for i in range(len(trial["cases"]))]
    result = {
        "kind": "global_response_upper_audit_v1",
        "source": str(SRC / "trial_discovery.json"),
        "full_sector_enumeration_avoided_by_this_audit": False,
        "reason": "The response vectors are generated from inherited upper witnesses; replay is exact but discovery is not independent.",
        "rows": rows,
        "decision": "Do not claim a compact upper solver from these trials. The response-lifted witness remains a valid variational upper certificate, but its H6 result retains 200-state inherited support and its fresh-transfer result retains 54-state support from an HF seed.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
