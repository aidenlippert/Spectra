"""Measure two live development branches, then freeze a policy prior before campaigns."""
import json
from pathlib import Path
import sys

from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.v2.campaign import generate, evaluate, score, read, write, model_cost
from research.rsi_discovery_20260916.v2.discovery import prompt, primitive_manifest
from research.rsi_discovery_20260916.v2.algebra import BASE_SOURCE


if __name__ == "__main__":
    acquired, pilot, output = map(Path, sys.argv[1:])
    primitive, parent = read(acquired / "proposal.json"), read(pilot / "proposal.json")
    primitive_check, parent_check = read(acquired / "validation.json"), read(pilot / "validation.json")
    tasks = [r["task"] for r in parent_check["tasks"]]
    ledger = Ledger(output)
    # Separate easier initial-state training from the preserved failed advanced-parent branches.
    parent = {"source": BASE_SOURCE}
    parent_check = evaluate(ledger, {"kind": "map", "source": BASE_SOURCE, "tasks": tasks, "repeats": 3}, "reference_parent")
    actions = ("continue_search", "combine_methods")
    jobs = [(output / action, prompt(tasks, [primitive_manifest(primitive, primitive_check)], action,
        parent["source"], [parent_check])) for action in actions]
    generate(jobs, parallel=2)
    archive = []
    for action in actions:
        folder = output / action
        proposal = read(folder / "proposal.json") if (folder / "proposal.json").exists() else None
        result = evaluate(ledger, {"kind": "map", "source": proposal["source"], "primitive": primitive["source"],
            "tasks": tasks, "repeats": 3}, action) if proposal else {"status": "implementation_failure"}
        costs = model_cost(folder)
        accepted = result["status"] == "verified_encoded_claim"
        archive.append({"before": "verified_encoded_claim", "action": action, "measured": True,
            "accepted": accepted, "cost_seconds": costs["metered_wall_seconds"] + result.get("complete_worker_seconds", 0),
            "construction_seconds": score(result) if accepted else None,
            "improves_parent": accepted and score(result) < score(parent_check),
            "parent_source_sha256": digest(parent["source"]), "source_sha256": digest(proposal["source"]) if proposal else None,
            "claim": "observed development branch only; future branches require live evaluation"})
    write(output / "archive.json", archive)
    write(output / "result.json", {"archive": archive, "audit": ledger.audit(), "costs": ledger.costs()})
    print(json.dumps(archive), flush=True)
