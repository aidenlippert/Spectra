"""One predeclared follow-up after the negative first RSI comparison.

Six fresh discovery campaigns are conditionally independent given the fixed
historical library. They are not six independent acquisitions of that library.
The first study's failures and complete acquisition costs remain explicit.
"""
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import time

from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.process import ROOT, run
from research.rsi_discovery_20260916.v2.campaign import generate, model_cost, write, read, score, ARMS
from research.rsi_discovery_20260916.v2.compound import prompt, BASE_SOURCE
from research.rsi_discovery_20260916.v2.library import MethodLibrary
from research.rsi_discovery_20260916.v2.records import Research, Method, State, Action, environment, choose_action

HERE = Path(__file__).resolve().parent


def tasks(index, transfer=False):
    seed = (702193 if transfer else 18317) + 1031 * index
    if not transfer:
        return [{"seed": seed, "base": {"family": "collective", "spatial": 5, "weight": -1},
                 "unique_words": 48, "operators": 20, "terms": 5},
                {"seed": seed + 19, "base": {"family": "chain", "spatial": 6, "width": 3, "weight": -1},
                 "unique_words": 64, "operators": 24, "terms": 8, "relations": True}]
    return [{"seed": seed, "base": {"family": "star", "spatial": 8, "weight": -1},
             "unique_words": 96, "operators": 40, "terms": 10, "relations": True},
            {"seed": seed + 29, "base": {"family": "chain", "spatial": 10, "width": 4, "weight": -1},
             "unique_words": 120, "operators": 48, "terms": 12, "zero_operator": True},
            {"seed": seed + 47, "base": {"family": "disjoint", "spatial": 12, "weight": -1,
             "mode_offset": 59, "mode_stride": 3}, "unique_words": 80, "operators": 32,
             "terms": 8, "large_integers": True, "relations": True},
            {"seed": seed + 61, "base": {"family": "collective", "spatial": 10, "weight": -1, "adjoint": True},
             "unique_words": 100, "operators": 40, "terms": 4, "zero_operator": True}]


def evaluate(ledger, request, label, seconds=180):
    folder = ledger.root / "evaluations" / label
    folder.mkdir(parents=True)
    write(folder / "request.json", request)
    with ledger.measure("compound_exact_evaluation", label_id=label, request_sha256=digest(request)) as cost:
        process = run([sys.executable, "-B", "-S", "-m", "research.rsi_discovery_20260916.v2.compound",
            str(folder / "request.json"), str(folder / "result.json")], folder / "process", seconds)
        cost.update(process)
        cost["stdout_sha256"] = ledger.blob(Path(process["stdout"]).read_bytes())
        cost["stderr_sha256"] = ledger.blob(Path(process["stderr"]).read_bytes())
        if process["exit_code"] or not (folder / "result.json").exists():
            cost["status"] = "timeout" if process["timeout"] else "failed"
            return {"status": State.COST.value if process["timeout"] else State.REPAIRABLE.value, "process": process}
        result = read(folder / "result.json")
        cost["scientific_state"] = result["status"]
        cost["output_sha256"] = ledger.blob((folder / "result.json").read_bytes())
    return result


def unavailable_state(folder):
    ledger = Ledger(folder)
    rows = ledger.events("operation_finished")
    timeout = any(row["result"].get("timeout") for row in rows)
    return {"status": State.COST.value if timeout else State.REPAIRABLE.value,
            "reason": "No proposal returned within the fixed native model budget" if timeout else "Native model response unavailable"}


def main(previous, conformance, output):
    output.mkdir(parents=True, exist_ok=False)
    start = time.time()
    before = read(previous / "result.json")
    if before["predeclared_operational_test_passed"]:
        raise ValueError("This follow-up protocol was specified after a negative first test")
    qualified = {(r["campaign"], r["arm"]) for r in read(conformance / "result.json")["results"]
                 if r["status"] == State.ACCEPTED.value}
    candidates = []
    for path in previous.glob("campaign_*/result.json"):
        old = read(path)
        for arm, record in old["arms"].items():
            if record.get("method_id") and (path.parent.name, arm) in qualified:
                candidates.append((score(record["transfer"]), path.parent, arm, record))
    _, source_ledger, source_arm, source_arm_record = min(candidates, key=lambda row: row[0])
    library = MethodLibrary(Ledger(source_ledger))
    retained_id = source_arm_record["method_id"]
    retained = library.methods[retained_id]
    retained_ref = {"ledger": str(source_ledger), "method_id": retained_id}
    source = library.ledger.read_blob(retained["source_sha256"]).decode()
    # This chosen primitive is closed pure arithmetic; a dependent primitive needs its full closure imported.
    if "car_kernel(" in source:
        raise ValueError("Import the entire dependency closure before using this primitive")
    prior = read(previous.parent / "policy_archive.json")
    source_names = ("compound_campaign.py", "compound.py", "campaign.py", "model.py", "model_job.py", "records.py", "library.py", "algebra.py", "language.py")
    sources = {str((HERE / name).relative_to(ROOT)): hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in source_names}
    protocol = {"experiment": "polynomial_operator_map_followup", "campaigns": 6, "arms": ARMS,
        "model": "gpt-6-astra", "effort": "high", "proposal_calls_per_arm": 2, "proposal_timeout_seconds": 180,
        "evaluation_timeout_seconds": 180, "repetitions": 3, "proposal_parallelism": 4,
        "independence": "Six fresh discovery campaigns conditional on one shared, already acquired library; not independent library acquisitions.",
        "historical_acquisition": {"study": str(previous), "study_passed": False, "selection": "fastest exact transfer method that also passed all real operator-channel conformance",
            "campaign": source_ledger.name, "arm": source_arm, "method_id": retained_id, "source_sha256": retained["source_sha256"],
            "cost_charge": "The complete earlier study and subsequent conformance are acquisition costs, not just the winning proposal."},
        "new_problem": "Exact symmetric Gram maps of sparse polynomial operators, not monomial dictionaries",
        "library_exposure": "Operational contract and measured regime; full implementation not injected into prompts. The live callable is installed automatically.",
        "fixed_policy": "continue search", "replay_policy": "frozen previous measured action prior, with typed response to actual new observations",
        "selection": "accepted development construction time only; freeze before new transfer; no transfer repair",
        "primary": {"contrast": "D versus A", "all_six_wins": True, "median_speedup_min": 1.25,
                    "B_over_A_and_D_over_C_median_min": 1.2, "all_D_generated_and_exact_transfer": True,
                    "one_sided_sign_p_max": 0.025, "multiplicity": "Bonferroni family of two RSI studies; report both, including the negative first study"},
        "universal_RSI_claim": False, "cumulative_payback_requires_separate_full_accounting": True,
        "development": [tasks(i) for i in range(6)], "transfer": [tasks(i, True) for i in range(6)], "sources": sources}
    write(output / "protocol.json", protocol)
    contexts = []
    files = [HERE / "compound.py", HERE / "library.py", HERE / "algebra.py", HERE / "language.py",
             ROOT / "experiments/marginal_symbolic.py", ROOT / "experiments/marginal_hunt_car.py"]
    for index in range(6):
        ledger = Ledger(output / ("campaign_" + str(index + 1)))
        ledger.append("protocol_frozen", sha256=ledger.blob((output / "protocol.json").read_bytes()))
        env = environment({"development": tasks(index), "withheld_transfer_hash": digest(tasks(index, True))}, files, [])
        research = Research(ledger, {"claim": protocol["new_problem"], "campaign": index}, env)
        node = research.observe(action=Action.TRANSFER, state=State.ACCEPTED,
            claim="Existing monomial-map method and all 146 actual operator-channel checks reopened from immutable history",
            evidence={"origin": retained_ref, "source_sha256": retained["source_sha256"], "prior_result": str(previous / "result.json"),
                      "conformance": str(conformance / (source_ledger.name + "_" + source_arm + ".json"))})
        policy = research.observe(action=Action.SEARCH, state=State.ACCEPTED,
            claim="Historical policy observations only; new-domain outcomes require live computation", evidence={"archive": prior})
        primitive = Method("monomial_map", ledger.blob(source), {"words": "homogeneous-charge degree<=3 CAR monomials", "cache": "mutable owned namespace"},
            {"map": "canonical word -> upper Gram coordinate -> exact integer"}, tuple(retained["assumptions"]),
            ("monomial-map",), (node,), retained["cost_regime"], (), tuple(retained["failure_cases"]),
            "monomial_map({'words': words, 'cache': payload['cache']})", "enabling")
        method_id = research.retain(primitive)
        contexts.append({"index": index, "ledger": ledger, "research": research, "import": node, "policy": policy,
            "primitive_id": method_id, "arms": {a: {"attempts": [], "reads": [policy], "best": None, "score": math.inf,
                                                  "last_source": BASE_SOURCE} for a in ARMS}})
    for round_index in range(2):
        jobs, pending = [], []
        for c in contexts:
            for arm, flags in ARMS.items():
                a = c["arms"][arm]
                available = set(a["reads"] + ([c["import"]] if flags["acquired"] else []))
                methods = c["research"].applicable("monomial-map", available)
                manifest = [] if not methods else [{"name": "monomial_map", "call": primitive.usage,
                    "inputs": primitive.inputs, "outputs": primitive.outputs, "assumptions": primitive.assumptions,
                    "failure_cases": primitive.failure_cases, "cost_regime": "Cold monomial-map construction verified at 0.471 aggregate seconds on the previous four transfer dictionaries; no timing claim on this new task.",
                    "scope": "Full exact symmetrized monomial Gram map; returns all canonical rows and (i,j) coefficients, with no scalar/sector quotient."}]
                state = a["attempts"][-1]["validation"]["status"] if a["attempts"] else State.OPEN.value
                action = Action.OPEN if not round_index else choose_action(state, prior, fixed=not flags["replay"])
                parent = a["last_source"] if state in (State.REPAIRABLE.value, State.REFUTED.value) else a["best"]["source"] if a["best"] else BASE_SOURCE
                feedback = [{"prior_action_observations_from_another_domain": prior}] + [v["validation"] for v in a["attempts"]]
                folder = c["ledger"].root / ("arm_" + arm) / ("model_" + str(round_index + 1))
                jobs.append((folder, prompt(tasks(c["index"]), manifest, action.value, parent, feedback)))
                pending.append((c, arm, action, available, folder))
        print(json.dumps({"stage": "new_problem_discovery", "round": round_index + 1, "jobs": len(jobs)}), flush=True)
        generate(jobs, timeout=180)
        for c, arm, action, available, folder in pending:
            a = c["arms"][arm]
            proposed = read(folder / "proposal.json") if (folder / "proposal.json").exists() else None
            if proposed:
                a["last_source"] = proposed["source"]
                value = evaluate(c["ledger"], {"source": proposed["source"], "tasks": tasks(c["index"]), "repeats": 3,
                    "retained": retained_ref if ARMS[arm]["acquired"] else None}, arm + "_" + str(round_index))
            else:
                value = unavailable_state(folder)
            node = c["research"].observe(action=action, state=value["status"], reads=tuple(sorted(available)),
                claim="new polynomial-operator map algorithm", evidence=value, cost=model_cost(folder))
            a["reads"].append(node)
            a["attempts"].append({"source_sha256": c["ledger"].blob(proposed["source"]) if proposed else None,
                "validation": value, "node": node, "action": action.value, "cost": model_cost(folder)})
            if score(value) < a["score"]:
                a.update(best=proposed, score=score(value), best_node=node)
            print(json.dumps({"stage": "development_result", "campaign": c["index"] + 1, "arm": arm,
                "round": round_index + 1, "status": value["status"], "seconds": score(value) if math.isfinite(score(value)) else None}), flush=True)
    for c in contexts:
        for arm, a in c["arms"].items():
            source = a["best"]["source"] if a["best"] else BASE_SOURCE
            a["frozen_source"] = c["ledger"].blob(source)
            a["generated_method_found"] = a["best"] is not None and source != BASE_SOURCE
            dependencies = a["reads"] + ([c["import"]] if ARMS[arm]["acquired"] else [])
            if a["best"]:
                retained_new = Method("polynomial_map_" + arm, a["frozen_source"], {"operators": "sparse exact polynomial dictionary", "cache": "empty"},
                    {"map": "exact symmetric Gram coefficient map"}, (a["best"]["assumptions"], "fresh inputs require exact checking"),
                    ("polynomial-map",), (a["best_node"],), {"development_seconds": a["score"]}, tuple(dependencies),
                    (a["best"]["failure_cases"],), "propose({'operators': operators, 'cache': {}})", "enabling")
                a["method_id"] = c["research"].retain(retained_new)
            write(c["ledger"].root / ("arm_" + arm) / "frozen.json", {"source": source, "source_sha256": a["frozen_source"],
                  "dependencies": dependencies, "generated_method_found": a["generated_method_found"]})
            c["ledger"].append("method_frozen_before_transfer", arm=arm, source_sha256=a["frozen_source"], dependencies=dependencies)
    print(json.dumps({"stage": "frozen_compound_transfer"}), flush=True)
    for c in contexts:
        index = c["index"]
        order = list(ARMS)[index % 4:] + list(ARMS)[:index % 4]
        for arm in order:
            a = c["arms"][arm]
            env = environment({"transfer": tasks(index, True)}, files, [a.get("method_id"), retained_id if ARMS[arm]["acquired"] else None])
            c["research"].env = env["fingerprint"]
            c["ledger"].append("live_environment_transition", environment=env)
            value = evaluate(c["ledger"], {"source": c["ledger"].read_blob(a["frozen_source"]).decode(),
                "tasks": tasks(index, True), "repeats": 3, "retained": retained_ref if ARMS[arm]["acquired"] else None}, "transfer_" + arm)
            a["transfer"] = value
            c["research"].observe(action=Action.TRANSFER, state=value["status"], reads=a["reads"],
                claim="unchanged polynomial-map method on fresh representations and integer scales", evidence=value)
            print(json.dumps({"stage": "transfer_result", "campaign": index + 1, "arm": arm, "status": value["status"],
                              "seconds": score(value) if math.isfinite(score(value)) else None}), flush=True)
        write(c["ledger"].root / "result.json", {"index": index, "arms": c["arms"], "costs": c["ledger"].costs(), "audit": c["ledger"].audit()})
    comparisons = []
    for c in contexts:
        times = {a: score(v["transfer"]) for a, v in c["arms"].items()}
        valid = all(math.isfinite(v) and v > 0 for v in times.values())
        ratios = {x + "/" + y: times[x] / times[y] for x, y in (("A", "D"), ("A", "B"), ("C", "D"), ("A", "C"), ("B", "D"))} if valid else {}
        comparisons.append({"campaign": c["index"] + 1, "comparable": valid, "speedups": ratios,
                            "construction_seconds": {a: t if math.isfinite(t) else None for a, t in times.items()}})
    comparable = all(row["comparable"] for row in comparisons)
    medians = {key: statistics.median(row["speedups"][key] for row in comparisons) for key in comparisons[0]["speedups"]} if comparable else {}
    wins = sum(row.get("speedups", {}).get("A/D", 0) > 1 for row in comparisons)
    p = sum(math.comb(6, k) for k in range(wins, 7)) / 64 if comparable else None
    stable = all(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == key for name, key in sources.items())
    passed = comparable and stable and wins == 6 and p <= .025 and medians["A/D"] >= 1.25 and medians["A/B"] >= 1.2 and medians["C/D"] >= 1.2 and all(c["arms"]["D"]["generated_method_found"] for c in contexts)
    result = {"protocol_sha256": digest(protocol), "elapsed_wall_seconds": time.time() - start, "source_snapshot_unchanged": stable,
        "comparisons": comparisons, "median_speedups": medians, "D_over_A_wins": wins, "paired_sign_p_one_sided": p,
        "two_study_Bonferroni_p": min(1, 2 * p) if p is not None else None, "predeclared_operational_test_passed": passed,
        "earlier_study_passed": False, "universal_RSI_proved": False, "full_lifetime_payback_proved": False,
        "scope": "Conditional discovery improvement on a new representation problem using one prior acquired library; no indefinite self-improvement theorem"}
    write(output / "result.json", result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main(*(Path(x).resolve() for x in sys.argv[1:]))
