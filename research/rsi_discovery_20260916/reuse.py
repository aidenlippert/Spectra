"""Reopen acquired methods in a separate campaign with fresh exact checks."""
import json
from pathlib import Path
import time

from research.rsi_discovery_20260916 import candidates, lemmas
from research.rsi_discovery_20260916.campaign import attempt, generate_tasks, source_snapshot, check_sources, write, compare_arms
from research.rsi_discovery_20260916.ledger import Ledger, canonical


def load_library(directory):
    prior = Ledger(directory)
    audit = prior.audit()
    if audit["unfinished_operations"] or audit["unfinished_attempts"] or not prior.events("campaign_finished"):
        raise ValueError("Only completed, audited campaign libraries can be reused")
    finished = prior.events("campaign_finished")[-1]
    report = json.loads(prior.read_blob(finished["report_sha256"]))
    selected = next((r for r in prior.events("retained_method")
                     if r["category"] == "constructor" and r["name"] == "development_selected"), None)
    if selected is None:
        raise ValueError("No accepted development constructor available")
    attempts = prior.events("attempt_finished")
    if not any(a["status"] == "verified_instance" and a["source_sha256"] == selected["source_sha256"]
               and a["witness_sha256"] == selected["evidence"] for a in attempts):
        raise ValueError("Retained source has no accepted attempt evidence")
    proof = json.loads(prior.read_blob(selected["evidence"]))
    if proof["status"] != "verified_instance":
        raise ValueError("Unverified source evidence")
    # Reuse is only a proposal privilege. Every new scientific witness is checked
    # independently, regardless of this previous success.
    kernel = next((r for r in prior.events("retained_method") if r["category"] == "kernel"), None)
    accelerated = bool(kernel and prior.read_blob(kernel["source_sha256"]).decode() == candidates.GRAM_FLINT)
    return {"constructor_source": prior.read_blob(selected["source_sha256"]).decode(),
            "constructor_sha256": selected["source_sha256"], "accelerated_kernel": accelerated,
            "prior_report_sha256": finished["report_sha256"], "prior_wall_seconds": report["total_run_wall_seconds"],
            "prior_cost_exclusions": report["known_cost_exclusions"]}


def reuse_accounting(rows, prior_seconds, elapsed, generation_seconds):
    # Shared fresh input generation cancels between the two arms. Charge all
    # remaining non-attempt work (library/proof replay and controller overhead)
    # to acquisition, rather than silently dropping it from the comparison.
    reactivation = max(0., elapsed - sum(r["cost_seconds"] for r in rows) - generation_seconds)
    return {**compare_arms(rows, "without_acquired_kernel", "with_acquired_kernel", prior_seconds + reactivation),
            "prior_run_seconds": prior_seconds, "library_reactivation_and_overhead_seconds": reactivation,
            "shared_generation_seconds": generation_seconds}


def run_saved(prior_path, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    ledger = Ledger(output)
    started = time.monotonic()
    with ledger.measure("method_library_load") as metered:
        methods = load_library(prior_path)
        metered.update(prior_report_sha256=methods["prior_report_sha256"])
    sources = source_snapshot(ledger)
    ledger.append("inherited_acquisition", **methods)
    # Fresh geometries declared before loading any candidate outcomes.
    specs = [
        {"name": "h4_asymmetric_future", "split": "validation", "basis": "sto-3g",
         "geometry": [["H", [0., 0., z]] for z in [0., 1.0, 2.25, 3.55]]},
        {"name": "water_future", "split": "validation", "basis": "sto-3g",
         "geometry": [["O", [0., 0., 0.]], ["H", [.72, .04, .63]], ["H", [-.79, -.02, .59]]]},
    ]
    tasks = generate_tasks(ledger, sources, specs)
    # Recompile the exact templates on this environment; old proof receipts are
    # evidence, not a substitute for checking a newly enabled proof capability.
    lemmas.discover(ledger)
    enabled = methods["accelerated_kernel"] and any(r.get("name") == "cross_term_bounds" and r["status"] == "general_lemma"
                                                   for r in ledger.events("lemma"))
    rows = []
    for index, task in enumerate(tasks):
        options = [("current_procedure", candidates.constructor(), False),
                   ("retained_constructor", methods["constructor_source"], False),
                   ("without_acquired_kernel", candidates.adaptive_constructor(False), False),
                   ("with_acquired_kernel", candidates.adaptive_constructor(enabled), enabled)]
        if index % 2:
            options.reverse()
        for name, source, acquired in options:
            rows.append(attempt(ledger, task, source, name, acquired=acquired, role="separate_future_campaign"))
    check_sources(sources)
    total = {name: sum(r["cost_seconds"] for r in rows if r["family"] == name) for name, _, _ in options}
    elapsed = time.monotonic() - started
    generation_seconds = sum(r["wall_seconds"] for r in ledger.events("operation_finished") if r["label"] == "molecular_input_generation")
    report = {"source_run": str(Path(prior_path).resolve()), "source_report_sha256": methods["prior_report_sha256"],
              "attempts": rows, "totals": total, "prior_acquisition_seconds": methods["prior_wall_seconds"],
              "kernel_ablation": reuse_accounting(rows, methods["prior_wall_seconds"], elapsed, generation_seconds),
              "scientific_speedup_established": False, "costs": ledger.costs(),
              "known_cost_exclusions": methods["prior_cost_exclusions"],
              "total_run_wall_seconds": elapsed, "audit": ledger.audit()}
    write(output / "report.json", report)
    ledger.append("campaign_finished", report_sha256=ledger.blob(canonical(report)))
    print(json.dumps({"report": str(output / "report.json"), "totals": total,
                      "kernel_ablation": report["kernel_ablation"]}), flush=True)
    return report
