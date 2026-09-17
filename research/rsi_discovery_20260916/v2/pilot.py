import json
from pathlib import Path
import sys

from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.discovery import prompt, primitive_manifest
from research.rsi_discovery_20260916.v2.model import propose


def main():
    acquired, out = map(Path, sys.argv[1:])
    primitive = json.loads((acquired / "proposal.json").read_text())
    validation = json.loads((acquired / "validation.json").read_text())
    tasks = [{"family": "collective", "spatial": 6, "weight": -1},
             {"family": "chain", "spatial": 7, "width": 3, "weight": -1, "permutation_seed": 438}]
    ledger = Ledger(out)
    value = propose(ledger, prompt(tasks, [primitive_manifest(primitive, validation)], "open_mechanism"), "pilot-map", timeout=150)
    (out / "proposal.json").write_text(json.dumps(value, indent=2) + "\n")
    if value:
        (out / "request.json").write_text(json.dumps({"kind": "map", "source": value["source"],
            "primitive": primitive["source"], "tasks": tasks, "repeats": 3}))
    print(json.dumps({"proposal": bool(value), "costs": ledger.costs(), "audit": ledger.audit()}))


if __name__ == "__main__":
    main()
