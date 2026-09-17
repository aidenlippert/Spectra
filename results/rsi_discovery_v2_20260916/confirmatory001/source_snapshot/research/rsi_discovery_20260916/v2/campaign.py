"""Predeclared independent acquisition/discovery/transfer campaigns, four ablations.

Only proposal inference runs concurrently. Scientific evaluations run serially;
final transfer runs after every proposal has finished. Model contexts are isolated.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import time

from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.process import ROOT, run
from research.rsi_discovery_20260916.v2.acquire import PROMPT as ACQUIRE
from research.rsi_discovery_20260916.v2.algebra import BASE_SOURCE
from research.rsi_discovery_20260916.v2.discovery import prompt, primitive_manifest
from research.rsi_discovery_20260916.v2.records import (Research, Method, State, Action, environment, choose_action)

HERE = Path(__file__).resolve().parent
ARMS = {"A": {"acquired": False, "replay": False}, "B": {"acquired": True, "replay": False},
        "C": {"acquired": False, "replay": True}, "D": {"acquired": True, "replay": True}}


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    def finite(v):
        if isinstance(v, float) and not math.isfinite(v):
            return None
        if isinstance(v, dict):
            return {k: finite(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [finite(x) for x in v]
        return v
    path.write_text(json.dumps(finite(value), indent=2, allow_nan=False) + "\n")


def read(path):
    return json.loads(path.read_text())


def development(index):
    return [{"family": "collective", "spatial": 5 + index % 2, "weight": -1,
             "permutation_seed": 8801 + index * 79},
            {"family": "chain", "spatial": 6 + index % 2, "width": 3,
             "weight": -1, "adjoint": bool(index % 2), "permutation_seed": 127 + index * 93}]


def transfer(index):
    return [{"family": "star", "spatial": 7, "weight": -1, "permutation_seed": 99001 + index},
            {"family": "disjoint", "spatial": 15, "weight": -1, "adjoint": True,
             "mode_offset": 101, "mode_stride": 3, "permutation_seed": 773 + index},
            {"family": "collective", "spatial": 10, "weight": -1, "adjoint": True,
             "permutation_seed": 67931 + index},
            {"family": "chain", "spatial": 10, "width": 4, "weight": -3,
             "mode_offset": 71, "mode_stride": 2, "permutation_seed": 319 + index}]


def score(result):
    if result.get("status") != State.ACCEPTED.value:
        return math.inf
    return sum(row["median_seconds"] for row in result.get("tasks", []))


def evaluate(ledger, request, label, timeout=90):
    folder = ledger.root / "evaluations" / label
    folder.mkdir(parents=True)
    write(folder / "request.json", request)
    with ledger.measure("scientific_evaluation", label_id=label, request_sha256=digest(request)) as cost:
        result = run([sys.executable, "-B", "-m", "research.rsi_discovery_20260916.v2.evaluate",
                      str(folder / "request.json"), str(folder / "result.json")], folder / "process", timeout)
        cost.update(result)
        if result["exit_code"] or not (folder / "result.json").exists():
            cost["status"] = "timeout" if result["timeout"] else "failed"
            return {"status": State.COST.value if result["timeout"] else State.REPAIRABLE.value,
                    "process": result}
        value = read(folder / "result.json")
        cost["scientific_state"] = value["status"]
        cost["output_sha256"] = ledger.blob((folder / "result.json").read_bytes())
    return value


def generate(jobs, parallel=4):
    """No shared ledger/process CPU counters across concurrent model calls."""
    def one(job):
        directory, text = job
        request = directory.parent / (directory.name + "_request.json")
        write(request, {"prompt": text, "timeout": 150})
        log = directory.parent / (directory.name + "_driver.log")
        with log.open("wb") as stream:
            try:
                process = subprocess.run([sys.executable, "-B", "-m", "research.rsi_discovery_20260916.v2.model_job",
                                          str(request), str(directory)], cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                         timeout=175)
            except subprocess.TimeoutExpired:
                return str(directory), 124
        return str(directory), process.returncode
    with ThreadPoolExecutor(max_workers=parallel) as pool:
        # All results are inspected by the caller; failed jobs remain explicit observations.
        return list(pool.map(one, jobs))


def model_cost(folder):
    return read(folder / "receipt.json")["costs"] if (folder / "receipt.json").exists() else {
        "metered_wall_seconds": 0, "unmeasured_failure": True}


def run_campaign(output, count=5, policy_archive=()):
    output.mkdir(parents=True, exist_ok=False)
    started = time.time()
    source_names = ("campaign.py", "model.py", "model_job.py", "records.py", "algebra.py", "evaluate.py", "language.py", "discovery.py", "acquire.py")
    sources = {str((HERE / name).relative_to(ROOT)): hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in source_names}
    protocol = {"version": 2, "independent_campaigns": count, "arms": ARMS,
        "generation_one": "one independent model proposal per campaign, exact CAR acquisition gate",
        "generation_two": "two model proposals per arm; select only by accepted development construction time",
        "generation_three": "freeze source, four unseen support/size/adjoint/relabeling regimes; no repair from transfer",
        "model": "gpt-6-astra", "reasoning_effort": "high", "model_timeout_seconds": 150,
        "evaluation_timeout_seconds": 90, "repetitions": 3, "proposal_parallelism": 4,
        "execution_order": "all inference completed before serial transfer; arm order rotated by campaign",
        "target": "exact homogeneous-charge degree-three SOS Gram map including every cross term",
        "primary_contrast": "D versus A paired campaign-level transfer construction time",
        "primary_acceptance": {"all_D_transfer_exact": True, "all_five_paired_wins_required": True,
            "median_D_speedup_over_A_at_least": 1.2, "D_must_be_new_generated_method": True,
            "both_B_over_A_and_D_over_C_median_at_least": 1.2},
        "policy_contrast": "C/A and D/B reported separately; no policy-improvement claim unless supported",
        "statistics": "one-sided paired sign test on five campaign units; descriptive secondary contrasts, no candidate pseudo-replication",
        "costs": "acquisition, proposal, activation, construction, all validation and failures; marginal and cumulative separate",
        "cost_exclusions": ["parent implementation tokens unavailable", "setup prior to this protocol", "subscription dollar allocation unknown"],
        "universal_RSI_claim": False, "policy_training_archive": list(policy_archive),
        "development": [development(i) for i in range(count)], "transfer": [transfer(i) for i in range(count)],
        "sources": sources}
    write(output / "protocol.json", protocol)
    campaigns = []
    for index in range(count):
        ledger = Ledger(output / ("campaign_" + str(index + 1)))
        ledger.append("protocol_frozen", protocol_sha256=ledger.blob((output / "protocol.json").read_bytes()))
        evaluator_files = [HERE / "algebra.py", HERE / "evaluate.py", HERE / "language.py",
                           ROOT / "experiments/marginal_symbolic.py", ROOT / "experiments/marginal_hunt_car.py"]
        env = environment({"development": development(index), "transfer_hash": digest(transfer(index))}, evaluator_files, [])
        research = Research(ledger, {"claim": protocol["target"], "campaign": index}, env)
        policy_node = research.observe(action=Action.SEARCH, state=State.ACCEPTED,
            claim="frozen observations from separate policy-development branches; no unvisited outcomes",
            evidence={"archive": list(policy_archive)})
        campaigns.append({"ledger": ledger, "research": research, "arms": {}, "index": index, "policy_node": policy_node})
    jobs = [(c["ledger"].root / "acquisition", ACQUIRE) for c in campaigns]
    print(json.dumps({"stage": "independent_acquisitions", "jobs": len(jobs)}), flush=True)
    generate(jobs)
    for c in campaigns:
        ledger, research, index = c["ledger"], c["research"], c["index"]
        path = ledger.root / "acquisition/proposal.json"
        proposal = read(path) if path.exists() else None
        result = evaluate(ledger, {"kind": "kernel", "source": proposal["source"], "seed": 99197 + index}, "acquisition") if proposal else {"status": State.REPAIRABLE.value}
        node = research.observe(action=Action.OPEN, state=result["status"], claim="acquired CAR primitive on scoped exact gates",
                                evidence=result, cost=model_cost(ledger.root / "acquisition"))
        c.update(primitive=proposal if result["status"] == State.ACCEPTED.value else None,
                 acquisition_result=result, acquisition_node=node)
        if c["primitive"]:
            source = ledger.blob(proposal["source"])
            method = Method("car_kernel", source, {"payload": "left/right/cache"}, {"polynomial": "exact integral CAR"},
                (proposal["assumptions"], "universal candidate program correctness is not established"),
                ("CAR", "SOS-map"), (node,), {"measurements": result["benchmarks"]}, (),
                (proposal["failure_cases"],), "car_kernel({'left': l, 'right': r, 'cache': cache})", "enabling")
            c["primitive_id"] = research.retain(method)
        for arm in ARMS:
            c["arms"][arm] = {"attempts": [], "reads": [c["policy_node"]], "best": None, "best_score": math.inf}
    for round_index in range(2):
        jobs, pending = [], []
        for c in campaigns:
            for arm, flags in ARMS.items():
                a, r = c["arms"][arm], c["research"]
                available = set(a["reads"])
                if flags["acquired"] and c["primitive"]:
                    available.add(c["acquisition_node"])
                selected = r.applicable("SOS-map", available)
                methods = [primitive_manifest(c["primitive"], c["acquisition_result"])] if selected else []
                state = a["attempts"][-1]["validation"]["status"] if a["attempts"] else State.OPEN.value
                action = Action.OPEN if not round_index else choose_action(state, policy_archive, fixed=not flags["replay"])
                feedback = [{"common_policy_training_observations": list(policy_archive)}] + [v["validation"] for v in a["attempts"]]
                previous = a["best"]["source"] if a["best"] else None
                folder = c["ledger"].root / ("arm_" + arm) / ("model_" + str(round_index + 1))
                jobs.append((folder, prompt(development(c["index"]), methods, action.value, previous, feedback)))
                pending.append((c, arm, action, folder, available))
        print(json.dumps({"stage": "algorithm_discovery", "round": round_index + 1, "jobs": len(jobs)}), flush=True)
        generate(jobs)
        for c, arm, action, folder, available in pending:
            a, ledger = c["arms"][arm], c["ledger"]
            proposal = read(folder / "proposal.json") if (folder / "proposal.json").exists() else None
            primitive = c["primitive"]["source"] if ARMS[arm]["acquired"] and c["primitive"] else None
            request = {"kind": "map", "source": proposal["source"], "tasks": development(c["index"]),
                       "primitive": primitive, "repeats": 3} if proposal else None
            value = evaluate(ledger, request, arm + "_" + str(round_index)) if request else {"status": State.REPAIRABLE.value, "reason": "model proposal unavailable"}
            cost = model_cost(folder)
            node = c["research"].observe(action=action, state=value["status"], claim="generated direct symmetric map",
                reads=tuple(sorted(available)), evidence=value, cost=cost)
            a["reads"].append(node)
            item = {"source_sha256": ledger.blob(proposal["source"]) if proposal else None,
                    "validation": value, "cost": cost, "node": node, "action": action.value}
            a["attempts"].append(item)
            if score(value) < a["best_score"]:
                a.update(best=proposal, best_score=score(value), best_node=node)
            print(json.dumps({"stage": "development_result", "campaign": c["index"] + 1, "arm": arm,
                              "round": round_index + 1, "status": value["status"], "seconds": None if score(value) == math.inf else score(value)}), flush=True)
    # All selection is finished. Freeze every arm before any transfer response is visible.
    for c in campaigns:
        for arm, a in c["arms"].items():
            source = a["best"]["source"] if a["best"] else BASE_SOURCE
            a["frozen_source"] = c["ledger"].blob(source)
            a["generated_method_found"] = a["best"] is not None and source != BASE_SOURCE
            if a["best"]:
                dependencies = tuple(a["reads"] + ([c["acquisition_node"]] if ARMS[arm]["acquired"] and c["primitive"] else []))
                method = Method("symmetric_map_" + arm, a["frozen_source"],
                    {"payload.words": "arbitrary ordered homogeneous-charge CAR dictionary, degree<=3", "payload.cache": "empty dictionary"},
                    {"map": "canonical word -> upper Gram coordinate -> integer coefficient"},
                    (a["best"]["assumptions"], "new inputs require independent exact checking"),
                    ("SOS-map-constructor",), (a["best_node"],), {"development_seconds": a["best_score"]},
                    dependencies, (a["best"]["failure_cases"],), "propose({'words': words, 'cache': {}})", "enabling")
                a["method_id"] = c["research"].retain(method)
            write(c["ledger"].root / ("arm_" + arm) / "frozen.json", {"source": source,
                "source_sha256": a["frozen_source"], "generated_method_found": a["generated_method_found"], "reads": a["reads"]})
            c["ledger"].append("method_frozen_before_transfer", arm=arm, source_sha256=a["frozen_source"], dependencies=a["reads"])
    print(json.dumps({"stage": "frozen_transfer"}), flush=True)
    for c in campaigns:
        index = c["index"]
        order = list(ARMS)[index % 4:] + list(ARMS)[:index % 4]
        for arm in order:
            a = c["arms"][arm]
            primitive = c["primitive"]["source"] if ARMS[arm]["acquired"] and c["primitive"] else None
            source = c["ledger"].read_blob(a["frozen_source"]).decode()
            if a.get("method_id"):
                available = set(a["reads"] + ([c["acquisition_node"]] if ARMS[arm]["acquired"] and c["primitive"] else []))
                methods = c["research"].applicable("SOS-map-constructor", available)
                if a["method_id"] not in methods:
                    raise ValueError("Transfer attempted to use an unacquired method")
                source = c["ledger"].read_blob(methods[a["method_id"]].source_sha256).decode()
            env = environment({"transfer": transfer(index)}, evaluator_files,
                              [a.get("method_id"), c.get("primitive_id") if primitive else None])
            c["research"].env = env["fingerprint"]
            c["ledger"].append("live_environment_transition", environment=env,
                               reason="new inputs and retained method; old outcomes are not reused")
            value = evaluate(c["ledger"], {"kind": "map", "source": source, "primitive": primitive,
                "tasks": transfer(index), "repeats": 3}, "transfer_" + arm, 120)
            a["transfer"] = value
            c["research"].observe(action=Action.TRANSFER, state=value["status"], claim="unmodified method on withheld regimes",
                                   dependencies=a["reads"], evidence=value)
            print(json.dumps({"stage": "transfer_result", "campaign": index + 1, "arm": arm,
                "status": value["status"], "seconds": None if score(value) == math.inf else score(value)}), flush=True)
        write(c["ledger"].root / "result.json", {"index": index, "arms": c["arms"],
            "acquisition": c["acquisition_result"], "acquisition_cost": model_cost(c["ledger"].root / "acquisition"),
            "evaluation_costs": c["ledger"].costs(), "audit": c["ledger"].audit()})
    comparisons = []
    for c in campaigns:
        times = {a: score(v["transfer"]) for a, v in c["arms"].items()}
        valid = all(math.isfinite(t) and t > 0 for t in times.values())
        ratios = {x + "/" + y: times[x] / times[y] for x, y in (("A", "D"), ("A", "B"), ("C", "D"), ("A", "C"), ("B", "D"))} if valid else {}
        comparisons.append({"campaign": c["index"] + 1, "comparable": valid,
                            "construction_seconds": {a: t if math.isfinite(t) else None for a, t in times.items()}, "speedups": ratios})
    complete = all(r["comparable"] for r in comparisons)
    wins = sum(r.get("speedups", {}).get("A/D", 0) > 1 for r in comparisons)
    medians = {key: statistics.median(r["speedups"][key] for r in comparisons) for key in comparisons[0]["speedups"]} if complete else {}
    sign_p = sum(math.comb(count, k) for k in range(wins, count + 1)) / 2**count if complete else None
    operational = complete and count == 5 and wins == 5 and medians["A/D"] >= 1.2 and medians["A/B"] >= 1.2 and medians["C/D"] >= 1.2 and all(c["arms"]["D"]["generated_method_found"] for c in campaigns)
    unchanged = all(hashlib.sha256((ROOT / p).read_bytes()).hexdigest() == h for p, h in sources.items())
    result = {"protocol_sha256": digest(protocol), "elapsed_wall_seconds": time.time() - started,
        "source_snapshot_unchanged": unchanged, "comparisons": comparisons, "median_speedups": medians,
        "paired_sign_test_one_sided_p": sign_p, "predeclared_operational_test_passed": bool(operational and unchanged),
        "universal_RSI_proved": False, "cumulative_payback_established": False,
        "interpretation": "Scoped causal evidence requires the predeclared contrasts, exact transfer, complete costs, and independence limitations; kernel speed alone is insufficient."}
    write(output / "result.json", result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("output", type=Path)
    p.add_argument("--count", type=int, default=5)
    p.add_argument("--policy-archive", type=Path)
    a = p.parse_args()
    run_campaign(a.output.resolve(), a.count, read(a.policy_archive) if a.policy_archive else ())
