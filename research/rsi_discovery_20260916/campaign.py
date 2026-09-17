"""Executable discovery -> exact gates -> retention -> replay -> fresh reuse."""
import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import random
import sys
import time

from research.rsi_discovery_20260916 import candidates, lemmas, provider
from research.rsi_discovery_20260916.evaluator import constructor_payload
from research.rsi_discovery_20260916.ledger import Ledger, canonical, digest
from research.rsi_discovery_20260916.policy import Policy, train
from research.rsi_discovery_20260916.process import ROOT, run

PACKAGE = "research.rsi_discovery_20260916"
HERE = Path(__file__).resolve().parent


def write(path, value):
    Path(path).write_bytes(canonical(value) + b"\n")


def source_snapshot(ledger):
    dependencies = json.loads((HERE / "dependency_snapshot.json").read_text())
    sources = {}
    for record in dependencies["files"]:
        raw = (ROOT / record["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != record["sha256"]:
            raise ValueError("Inherited dependency changed: " + record["path"])
        sources[record["path"]] = ledger.blob(raw)
    for path in sorted(HERE.glob("*.py")):
        sources[str(path.relative_to(ROOT))] = ledger.blob(path.read_bytes())
    return sources


def check_sources(sources):
    for path, key in sources.items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != key:
            raise ValueError("Frozen source changed during campaign: " + path)


def specifications():
    def chain(count, spacing):
        return [["H", [0.0, 0.0, round(i * spacing, 8)]] for i in range(count)]
    return [
        {"name": "h2_development", "split": "development", "geometry": chain(2, 1.4), "basis": "sto-3g"},
        {"name": "h4_development", "split": "development", "geometry": chain(4, 1.4), "basis": "sto-3g"},
        {"name": "h4_changed_geometry", "split": "validation", "geometry": chain(4, 1.1), "basis": "sto-3g"},
        {"name": "water_nonchain", "split": "validation", "basis": "sto-3g",
         "geometry": [["O", [0., 0., 0.]], ["H", [.77, .03, .57]], ["H", [-.84, .01, .48]]]},
    ]


def generate_tasks(ledger, sources, specs=None):
    protocol = {"version": 1, "specifications": specs or specifications(), "target_width_Ha": "1/625",
                "sources": sources, "source_digest": digest(sources), "max_sector_labels": 2000,
                "attempt_timeout_seconds": 60, "serial_workers": 1,
                "validation_policy": "No selection updates from validation outcomes; report any reruns after implementation fixes",
                "kernel_repetitions": 3, "model_call_budget": 1,
                "cost_scope": "All automated proposal, generation, evaluation, failed attempts and policy development stages",
                "cost_exclusions": ["Earlier side-task setup", "Manual implementation and its model tokens", "Tests and proof preflight outside this run", "Parent-session tool overhead"],
                "money": "Existing subscription and local compute only"}
    write(ledger.root / "protocol.json", protocol)
    ledger.append("protocol_frozen", protocol_sha256=ledger.blob(canonical(protocol)))
    tasks = []
    for spec in protocol["specifications"]:
        directory = ledger.root / spec["name"]
        directory.mkdir()
        write(directory / "specification.json", spec)
        with ledger.measure("molecular_input_generation", task=spec["name"]) as cost:
            result = run([sys.executable, "-B", "-m", "research.transfer_solver_20260915.generate", str(directory)],
                         directory / "generate_process", 90)
            cost.update(result)
            if result["exit_code"]:
                raise RuntimeError("Molecular generation failed: " + result["stderr"])
        fixture = json.loads((directory / "fixture.json").read_text())
        fixture_hash = ledger.blob(canonical(fixture))
        task_id = digest({"specification": spec, "fixture": fixture_hash, "evaluator": digest(sources), "target": "1/625"})
        task = {"id": task_id, "name": spec["name"], "split": spec["split"], "fixture_sha256": fixture_hash,
                "fixture": fixture, "target": "1/625", "evaluator_sha256": digest(sources)}
        ledger.append("task_frozen", **{k: v for k, v in task.items() if k != "fixture"})
        tasks.append(task)
    return tasks


def attempt(ledger, task, source, family, parent=None, depth=0, payloads=None, repeats=1,
            proposal_seconds=0, acquired=False, role="development", timeout=60):
    started = time.monotonic()
    source_hash = ledger.blob(source)
    number = len(ledger.events("attempt_started"))
    node_id = f"{number:04d}-{source_hash[:12]}"
    directory = ledger.root / "attempts" / node_id
    directory.mkdir(parents=True)
    source_path = directory / "candidate.py"
    source_path.write_text(source)
    kernel = payloads is not None
    if not kernel:
        payload, allocation = constructor_payload(task["fixture"])
        payloads = [payload]
    else:
        allocation = {"factor_entries": sum(sum(map(len, p["factor"])) for p in payloads)}
    request = {"payloads": payloads, "repeats": repeats, "acquired_methods": acquired}
    write(directory / "request.json", request)
    ledger.append("attempt_started", id=node_id, task_id=task["id"], task_name=task["name"],
                  source_sha256=source_hash, parent=parent, depth=depth, family=family, role=role,
                  payload_sha256=ledger.blob(canonical(request)), allocation=allocation,
                  assumptions=["Fixed exact task", "60-second per-stage cap", "one local worker"])
    record = {"id": node_id, "task_id": task["id"], "task_name": task["name"], "parent": parent, "depth": depth,
              "family": family, "source_sha256": source_hash, "role": role,
              "status": "inconclusive", "utility": 0.0, "proposal_seconds": proposal_seconds}
    with ledger.measure("candidate_execution", attempt=node_id) as cost:
        execution = run([sys.executable, "-B", "-m", PACKAGE + ".worker", str(source_path),
                         str(directory / "request.json"), str(directory / "output.json")],
                        directory / "execution", timeout)
        cost.update(execution)
        if execution["exit_code"]:
            cost["status"] = "failed"
    record["execution"] = execution
    if not execution["exit_code"] and (directory / "output.json").exists():
        try:
            proposed = json.loads((directory / "output.json").read_text())
            record["candidate_output_sha256"] = ledger.blob(canonical(proposed))
            if kernel:
                check_input = {"payloads": payloads, "outputs": proposed["outputs"]}
            else:
                cert = proposed["outputs"][0]
                cert.update(kind="discovered_full_N_certificate_v1", fixture_sha256=digest(task["fixture"]))
                record["certificate_sha256"] = ledger.blob(canonical(cert))
                check_input = {"fixture": task["fixture"], "certificate": cert}
            write(directory / "check_input.json", check_input)
            with ledger.measure("independent_check", attempt=node_id) as cost:
                check = run([sys.executable, "-B", "-S", "-m", PACKAGE + ".evaluator", "grams" if kernel else "certificate",
                             str(directory / "check_input.json"), str(directory / "receipt.json")],
                            directory / "check", timeout)
                cost.update(check)
                if check["exit_code"]:
                    cost["status"] = "failed"
            record["check"] = check
            if not check["exit_code"]:
                receipt = json.loads((directory / "receipt.json").read_text())
                record.update(status="verified_instance", witness_sha256=ledger.blob(canonical(receipt)), receipt=receipt,
                              median_component_seconds=proposed["median_seconds"], samples_seconds=proposed["samples_seconds"])
                # Every accepted molecular result meets the same accuracy target.
                # Utility is a time objective; file size and numerical energy are not rewards.
                timing = proposed["median_seconds"] if kernel else execution["wall_seconds"] + check["wall_seconds"]
                record["utility"] = 1 / max(timing, 1e-9)
            else:
                record.update(status="rejected", failure_reason=Path(check["stderr"]).read_text()[-4000:])
        except (ValueError, KeyError, TypeError, IndexError, AttributeError, OverflowError) as error:
            record.update(status="rejected", failure_reason=str(error))
    else:
        record["failure_reason"] = Path(execution["stderr"]).read_text()[-4000:]
    record["cost_seconds"] = time.monotonic() - started + proposal_seconds
    record["directory"] = str(directory)
    ledger.append("attempt_finished", **record)
    print(json.dumps({"attempt": node_id, "task": task["name"], "family": family, "role": role,
                      "status": record["status"], "cost_seconds": round(record["cost_seconds"], 3)}), flush=True)
    return record


def kernel_payloads(cert, seed, validation=False):
    rng = random.Random(seed)
    factors = [r["factor"] for r in cert["blocks"] if len(r["factor"]) > 1]
    # Real certificate factors plus zero/sign/large-integer stress inputs.
    sizes = (1, 7, 32, 96) if not validation else (2, 11, 47)
    for n in sizes:
        factors.append([[rng.randint(-10**20, 10**20) if (i + j) % 3 else 0 for j in range(i + 1)] for i in range(n)])
    return [{"factor": f} for f in factors]


def explore_fresh(ledger, task, policy, label, budget_seconds=20):
    """Execute the same bounded proposal tree for each frozen policy arm."""
    options = [
        {"id": "0", "parent": None, "depth": 0, "family": "reference", "source": candidates.constructor()},
        {"id": "1", "parent": "0", "depth": 1, "family": "low_precision", "source": candidates.constructor(digits=3)},
        {"id": "2", "parent": "0", "depth": 1, "family": "subset", "source": candidates.constructor(subset=True)},
        {"id": "3", "parent": "1", "depth": 2, "family": "refinement",
         "source": candidates.constructor(digits=10, subset=True, guard="1/1000000")},
    ]
    completed, observed, real_ids = set(), [], {}
    spent, best, stale = 0., 0., 0
    while spent < budget_seconds and len(observed) < policy.max_attempts and stale < policy.patience:
        frontier = [{k: item[k] for k in ("id", "parent", "depth", "family")} for item in options
                    if item["id"] not in completed and (item["parent"] is None or item["parent"] in completed)]
        if not frontier:
            break
        for action in policy.choose(frontier, observed):
            if spent >= budget_seconds or len(observed) >= policy.max_attempts or stale >= policy.patience:
                break
            item = options[int(action["id"])]
            remaining = max(1, int(budget_seconds - spent))
            row = attempt(ledger, task, item["source"], item["family"],
                          parent=real_ids.get(item["parent"]), depth=item["depth"], role="fresh_policy_" + label,
                          timeout=min(60, remaining))
            observed.append(row)
            completed.add(item["id"])
            real_ids[item["id"]] = row["id"]
            spent += row["cost_seconds"]
            if row["utility"] > best:
                best, stale = row["utility"], 0
            else:
                stale += 1
    ledger.append("fresh_policy_result", task_id=task["id"], arm=label, policy=asdict(policy),
                  charged_seconds=spent, best_utility=best, attempts=[n["id"] for n in observed],
                  budget_seconds=budget_seconds, overshoot_seconds=max(0, spent - budget_seconds))
    return observed


def compare_arms(rows, left, right, acquisition):
    a, b = ([r for r in rows if r["family"] == arm] for arm in (left, right))
    matched = bool(a and b) and sorted(r["task_id"] for r in a) == sorted(r["task_id"] for r in b)
    complete = matched and all(r["status"] == "verified_instance" for r in a + b)
    left_cost, right_cost = sum(r["cost_seconds"] for r in a), sum(r["cost_seconds"] for r in b)
    return {"both_arms_accepted_on_same_tasks": complete, "left_attempted_seconds": left_cost,
            "right_attempted_seconds": right_cost, "acquisition_seconds": acquisition,
            "paid_cost_advantage_seconds": left_cost - right_cost - acquisition if complete else None,
            "comparison_status": "matched_finite_runs" if complete else "incomparable_failed_or_missing_arm"}


def run_campaign(output, use_codex=False):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    ledger = Ledger(output)
    wall_start = time.monotonic()
    sources = source_snapshot(ledger)
    versions = {k: importlib.metadata.version(k) for k in ("numpy", "scipy", "python-flint", "z3-solver", "pyscf")}
    ledger.append("environment", python=sys.version, platform=platform.platform(), packages=versions,
                  external_compute_purchased=False, prior_setup_costs_known=False)
    tasks = generate_tasks(ledger, sources)
    print("Frozen four molecular tasks and split before proposals.", flush=True)
    lemmas.discover(ledger)
    general = {r["name"]: r for r in ledger.events("retained_method") if r["category"] == "lemma"}
    dev = [t for t in tasks if t["split"] == "development"]
    histories, rows = [], []
    for task in dev:
        baseline = attempt(ledger, task, candidates.constructor(), "reference")
        bad = attempt(ledger, task, candidates.constructor(digits=3), "low_precision", parent=baseline["id"], depth=1)
        subset = attempt(ledger, task, candidates.constructor(subset=True), "subset", parent=baseline["id"], depth=1)
        refined = attempt(ledger, task, candidates.constructor(digits=10, subset=True, guard="1/1000000"),
                          "refinement", parent=bad["id"], depth=2)
        history = [baseline, bad, subset, refined]
        if task is dev[-1] and use_codex:
            proposal = provider.propose(ledger, candidates.constructor(),
                                        [{k: n.get(k) for k in ("family", "status", "cost_seconds", "failure_reason")} for n in history])
            if proposal:
                history.append(attempt(ledger, task, proposal["source"], "model_proposal", parent=baseline["id"], depth=1,
                                       proposal_seconds=proposal["cost_seconds"]))
        histories.append(history)
        rows += history
    dev_success = [r for r in histories[-1] if r["status"] == "verified_instance"]
    if not dev_success:
        raise RuntimeError("No accepted development constructor; no promotion allowed")
    fastest = min(dev_success, key=lambda r: r["execution"]["wall_seconds"] + r["check"]["wall_seconds"])
    ledger.append("retained_method", category="constructor", name="development_selected", source_sha256=fastest["source_sha256"],
                  evidence=fastest["witness_sha256"], status="benchmarked_algorithm", task_ids=[fastest["task_id"]],
                  scope="Accepted development model, selected by construction plus independent-check wall time")
    certificate = json.loads(ledger.read_blob(fastest["certificate_sha256"]))
    ktask = {"id": digest([dev[-1]["id"], "integer_gram", 917]), "name": "exact_gram_development"}
    kpayloads = kernel_payloads(certificate, 917)
    khistory = []
    for family, source in (("python_gram", candidates.GRAM_REFERENCE), ("triangular_gram", candidates.GRAM_TRIANGULAR),
                           ("flint_gram", candidates.GRAM_FLINT), ("invalid_gram", candidates.GRAM_BAD)):
        khistory.append(attempt(ledger, ktask, source, family, parent=khistory[0]["id"] if khistory else None,
                                depth=int(bool(khistory)), payloads=kpayloads, repeats=3))
    rows += khistory
    ksuccess = [r for r in khistory if r["status"] == "verified_instance"]
    kbest = min(ksuccess, key=lambda r: r["median_component_seconds"])
    use_acquired = kbest["family"] == "flint_gram" and "cross_term_bounds" in general
    ledger.append("retained_method", category="kernel", name=kbest["family"], evidence=kbest["witness_sha256"],
                  source_sha256=kbest["source_sha256"], status="benchmarked_algorithm", task_ids=[ktask["id"]],
                  scope="Finite development equivalence suite; always recheck new scientific certificates")
    # This real caller connects the learned exact component and checked arithmetic
    # lemma to a new constructor. The independent gate remains unchanged.
    adaptive = attempt(ledger, dev[-1], candidates.adaptive_constructor(use_acquired), "acquired_adaptive" if use_acquired else "adaptive",
                       parent=fastest["id"], depth=2, acquired=use_acquired)
    histories[-1].append(adaptive)
    rows.append(adaptive)
    policy, policy_records = train(ledger, histories, budget_seconds=20)
    training_costs = ledger.costs()
    check_sources(sources)
    ledger.append("development_frozen", policy=asdict(policy), constructor_source=fastest["source_sha256"],
                  acquisition_wall_seconds=time.monotonic() - wall_start, costs=training_costs,
                  acquired_margin_enabled=use_acquired)
    acquisition = time.monotonic() - wall_start
    future = []
    validation_kernel = []
    # Identical task/accuracy/resource envelopes; independently rerun every arm.
    # Alternate execution order by task to reduce simple warm-up order effects.
    for index, task in enumerate(t for t in tasks if t["split"] == "validation"):
        arms = [("current_procedure", candidates.constructor(), False),
                ("simple_adaptive", candidates.adaptive_constructor(False), False),
                ("with_acquired_methods", candidates.adaptive_constructor(use_acquired), use_acquired)]
        if index % 2:
            arms.reverse()
        for arm, source, acquired in arms:
            r = attempt(ledger, task, source, arm, acquired=acquired, role="fresh_validation")
            future.append(r)
        # Each arm gets an independent fresh history and identical proposal tree.
        for label, selected in (("fixed", Policy()), ("adaptive", Policy(name="adaptive", adaptive=True)),
                                ("replay_selected", policy)):
            explore_fresh(ledger, task, selected, label)
        successful = [n for n in future if n["task_id"] == task["id"] and n["status"] == "verified_instance"]
        if successful:
            cert = json.loads(ledger.read_blob(successful[0]["certificate_sha256"]))
            payloads = kernel_payloads(cert, 982 + index, validation=True)
            task_kernel = {"id": digest([task["id"], "kernel_validation"]), "name": task["name"] + "_gram"}
            for family, source in (("python_gram", candidates.GRAM_TRIANGULAR), ("retained_gram", ledger.read_blob(kbest["source_sha256"]).decode())):
                validation_kernel.append(attempt(ledger, task_kernel, source, family, payloads=payloads, repeats=3, role="fresh_validation"))
    check_sources(sources)
    totals = {}
    for r in future:
        totals[r["family"]] = totals.get(r["family"], 0) + r["cost_seconds"]
    current = totals.get("current_procedure", 0)
    acquired = totals.get("with_acquired_methods", 0)
    reuse = {**compare_arms(future, "current_procedure", "with_acquired_methods", acquisition),
             "acquisition_seconds": acquisition, "future_without_acquired_seconds": current,
             "future_with_acquired_seconds": acquired, "with_acquisition_seconds": acquisition + acquired,
             "all_arms_same_target": True, "replications_per_arm_per_task": 1,
             "scientific_speedup_established": False,
             "reason": "Two finite transfer tasks and one arm run each; acquisition is fully charged within the automated run"}
    for r in future:
        if r["status"] == "verified_instance":
            ledger.append("retained_method", category="constructor", name=r["family"], source_sha256=r["source_sha256"],
                          evidence=r["witness_sha256"], status="benchmarked_algorithm", task_ids=[r["task_id"]],
                          scope="Accepted finite rational model and measured run; no extrapolated theorem or speedup")
    report = {"protocol": str(output / "protocol.json"), "source_digest": digest(sources),
              "development_attempts": rows, "future_attempts": future, "kernel_validation": validation_kernel,
              "selected_policy": asdict(policy), "policy_baselines": policy_records[:2], "recursion": reuse,
              "fresh_policy_results": ledger.events("fresh_policy_result"),
              "costs": ledger.costs(), "total_run_wall_seconds": time.monotonic() - wall_start,
              "audit": ledger.audit(), "known_cost_exclusions": json.loads((output / "protocol.json").read_text())["cost_exclusions"]}
    write(output / "report.json", report)
    ledger.append("campaign_finished", report_sha256=ledger.blob(canonical(report)))
    print(json.dumps({"report": str(output / "report.json"), "recursion": reuse, "audit": report["audit"]}), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New directory; existing run directories are never overwritten")
    parser.add_argument("--provider", choices=("local", "codex"), default="local")
    parser.add_argument("--reuse-from", type=Path, help="Reopen a completed method library and run a separate fresh campaign")
    a = parser.parse_args()
    if a.reuse_from:
        from research.rsi_discovery_20260916.reuse import run_saved
        run_saved(a.reuse_from, a.output)
    else:
        run_campaign(a.output, use_codex=a.provider == "codex")


if __name__ == "__main__":
    main()
