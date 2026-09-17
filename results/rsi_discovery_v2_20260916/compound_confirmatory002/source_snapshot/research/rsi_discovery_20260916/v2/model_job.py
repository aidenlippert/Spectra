"""One independent native model invocation, with its own process accounting."""
import json
from pathlib import Path
import sys

from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.v2.model import propose


if __name__ == "__main__":
    request_path, directory = map(Path, sys.argv[1:])
    request = json.loads(request_path.read_text())
    ledger = Ledger(directory)
    value = propose(ledger, request["prompt"], "proposal", timeout=request.get("timeout", 150))
    (directory / "proposal.json").write_text(json.dumps(value, indent=2) + "\n")
    (directory / "receipt.json").write_text(json.dumps({"costs": ledger.costs(), "audit": ledger.audit()}, indent=2) + "\n")
