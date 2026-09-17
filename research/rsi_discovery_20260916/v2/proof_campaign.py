"""Metered mathematical gates and consumption of the new checked operation."""
import json
from pathlib import Path
import re
import statistics
import sys
import time

from experiments.marginal_symbolic import word_product
from research.rsi_discovery_20260916.ledger import Ledger
from research.rsi_discovery_20260916.process import run
from research.rsi_discovery_20260916.v2.algebra import oracle, gram_reference, task_words
from research.rsi_discovery_20260916.v2.certified import CertifiedTemplates
from research.rsi_discovery_20260916.v2.language import compile_program


def worker(proposal):
    source = json.loads(proposal.read_text())["source"]
    fn = compile_program(source, {"oracle": oracle})
    records = []
    for spatial in (6, 8, 10):
        words = task_words({"family": "collective", "spatial": spatial, "weight": -1})
        reference = gram_reference(words)
        for kind in ("original_exact_CAR", "checked_template_operation"):
            samples, receipt = [], None
            for _ in range(3):
                word_product.cache_clear()
                checked = CertifiedTemplates(fn)
                start = time.perf_counter()
                value = gram_reference(words, oracle if kind == "original_exact_CAR" else checked)
                samples.append(time.perf_counter() - start)
                if value != reference:
                    raise ValueError("Checked operation changed a coefficient map")
                receipt = checked.receipt()
            records.append({"spatial": spatial, "dictionary_words": len(words), "operation": kind,
                            "cold_seconds": samples, "median_seconds": statistics.median(samples), "receipt": receipt})
    print(json.dumps({"status": "verified_encoded_claim", "records": records}))


def main():
    if sys.argv[1] == "worker":
        worker(Path(sys.argv[2])); return
    proposal, output = map(Path, sys.argv[1:])
    ledger = Ledger(output)
    lean_source = Path(__file__).with_name("Transport.lean")
    with ledger.measure("transport_and_error_lemmas") as cost:
        result = run([str(Path.home() / ".elan/bin/lean"), str(lean_source)], output / "lean", 45)
        cost.update(result)
        text = Path(result["stdout"]).read_text()
        groups = re.findall(r"depends on axioms: \[([^\]]*)\]", text)
        axioms = {a.strip() for group in groups for a in group.split(",") if a.strip()}
        accepted = result["exit_code"] == 0 and len(groups) == 5 and axioms <= {"propext", "Quot.sound", "Classical.choice"}
        cost.update(source_sha256=ledger.blob(lean_source.read_bytes()), stdout_sha256=ledger.blob(text), accepted=accepted)
        if not accepted:
            raise ValueError("Lean transport/allowance proof gate failed")
    with ledger.measure("new_checked_operation_actual_map_consumption") as cost:
        result = run([sys.executable, "-B", "-S", "-m", "research.rsi_discovery_20260916.v2.proof_campaign",
                      "worker", str(proposal)], output / "consumer", 120)
        cost.update(result)
        if result["exit_code"]:
            raise ValueError("Checked-operation consumer failed")
        value = json.loads(Path(result["stdout"]).read_text())
        cost["receipt_sha256"] = ledger.blob(Path(result["stdout"]).read_bytes())
    summary = {"formal_lemmas": 5, "axioms": sorted(axioms), "checked_operation": value,
               "whole_Python_program_formally_verified": False, "costs": ledger.costs(), "audit": ledger.audit()}
    (output / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
