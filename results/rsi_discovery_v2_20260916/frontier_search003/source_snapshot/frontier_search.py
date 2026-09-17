"""Executable teacher-free representation proposals, with a separate fixed checker.

Candidate programs output only observable polynomials. They never choose the
Hamiltonian, initial state, controls, uncertainty budgets, target, or acceptance rule.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import traceback

from experiments.marginal_symbolic import decode, encode
from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.process import run
from research.rsi_discovery_20260916.v2.campaign import read, write
from research.rsi_discovery_20260916.v2.frontier import (
    frozen_problem, validate_observables, construct, replay, cadd, cscale, drift, coordinates)
from research.rsi_discovery_20260916.v2.language import compile_program
from research.rsi_discovery_20260916.v2.model import propose
from research.rsi_discovery_20260916.v2.records import Research, Method, State, Action, environment, choose_action

MODULE = "research.rsi_discovery_20260916.v2.frontier_search"
BASE_SOURCE = 'def propose(payload):\n    return [payload["W"], payload["Y"]]\n'


def inputs(folder):
    return read(folder / "fixture.json"), read(folder / "state.json")


def candidate_observables(source, data, state, max_additional=14):
    if type(max_additional) is not int or not 1 <= max_additional <= 14:
        raise ValueError("Observable budget must be an integer between one and fourteen")
    h, d, w, y = frozen_problem(data, state)
    capabilities = {"rational": F, "drift": drift, "cadd": cadd, "cscale": cscale,
                    "coordinates": coordinates}
    candidate = compile_program(source, capabilities)
    payload = {"H": h, "D": d, "W": w, "Y": y, "modes": data["modes"],
               "particles": data["particles"], "max_additional_observables": max_additional, "cache": {}}
    # The program never receives the MPS, a trajectory, or a full-state teacher.
    result = validate_observables(candidate(payload), data["modes"])
    if len(result) > max_additional:
        raise ValueError("Proposed basis exceeds the declared observable budget")
    return [{"real": encode(a), "imaginary": encode(b)} for a, b in result]


def checked_basis(encoded, modes):
    return validate_observables([(decode(item["real"], modes, 4),
                                  decode(item["imaginary"], modes, 4)) for item in encoded], modes)


def worker(action, request_path, output):
    request = read(request_path)
    data, state = inputs(Path(request["inputs"]))
    try:
        if action == "construct":
            result = {"status": "proposal_constructed", "observables":
                      candidate_observables(request["source"], data, state, request.get("max_additional", 14))}
        else:
            basis = checked_basis(read(Path(request["observables"]))["observables"], data["modes"])
            result = construct(data, state, basis)
            if action == "replay":
                replay(data, state, read(Path(request["certificate"])), basis)
            if any(name in sys.modules for name in ("numpy", "scipy", "quimb", "cvxpy", "pyscf")):
                raise AssertionError("Numerical package entered independent exact acceptance")
        write(output, result)
    except Exception as error:
        message = str(error)
        status = (State.COST if "envelope exceeded" in message else
                  State.REPRESENTATION if "linearly dependent" in message else State.REPAIRABLE)
        write(output, {"status": status.value,
                       "error": type(error).__name__ + ": " + str(error), "traceback": traceback.format_exc()})


PROMPT = """Discover a compact observable representation for Spectra's frozen H8 controlled dynamics.
The four-observable I,D,W,Y baseline has a rigorous but useless final D enclosure [-2,2].
We need D<=-3/5 for the original four half-unit phases. Both control amplitudes remain <=1/2.
Your program only proposes additional Hermitian observables; the checker always prepends I and D.
Input payload provides the exact Gaussian-rational CAR polynomials H,D,W,Y as pairs (real,imaginary),
plus modes=16, particles=8, max_additional_observables=14 and an empty cache. No state or trajectory is supplied.
Each polynomial is a dictionary canonical word -> exact Fraction coefficient. Words are tuples
of (creation_flag, mode); creators precede annihilators, each group sorted ascending.
Return a list of 1..14 polynomial pairs, each Hermitian, particle neutral and of degree <=4,
with at most 4096 terms per real/imaginary part. The complete basis including I,D must be independent.
Capabilities: rational(n,d) makes an exact Fraction; drift(G,O) computes i[G,O];
cadd(A,B), cscale(A,s); coordinates(O,basis) returns (coordinate_list, residual_dictionary).
Residual keys are (real_or_imaginary_part, canonical_word), with exact Fraction values.
coordinates requires independent basis elements. Lists/dicts can be built directly.
You can propose weighted commutator bases, localized correlated operators, alternative exact
coordinates, or a different sparse closure. You cannot change the physical task or norm gate.
The separate checker projects every i[H,O], i[D,O], i[W,O] exactly into the proposed basis,
reconstructs every residual coefficient, bounds residual norm by its coefficient l1 norm,
and propagates true dynamical-defect and rational-integrator errors over all four phases.
The growth bound uses the infinity norm of the resulting real coordinate matrix; enormous
coordinates are counterproductive. Only exact initial moments are contracted from the original MPS.
Optimize the rigorously accepted interval, not trajectory fit or a negative matrix eigenvalue.
The process has 120 seconds for basis construction and 180 seconds for exact acceptance.
Unresolved but correct defect bounds are retained as conditional research, never as a solved target.
Return one concise executable representation algorithm.
"""


def campaign(source_inputs, output, supplied=None, max_additional=14):
    output.mkdir(parents=True, exist_ok=False)
    ledger = Ledger(output)
    copied = output / "inputs"
    copied.mkdir()
    for name in ("fixture.json", "state.json"):
        (copied / name).write_bytes((source_inputs / name).read_bytes())
    data, state = inputs(copied)
    frozen_problem(data, state)
    files = [Path(__file__), Path(__file__).with_name("frontier.py"), Path(__file__).with_name("language.py")]
    env = environment({"fixture": digest(data), "state": digest(state), "max_additional_observables": max_additional}, files, [])
    research = Research(ledger, {"claim": "Original robust H8 task: D <= -3/5", "amplitude_limit": "1/2"}, env)
    protocol = {"model": "gpt-6-astra", "effort": "high", "proposal_calls": 0 if supplied else 2,
        "model_deadline_seconds": 150, "construction_deadline_seconds": 120, "checking_deadline_seconds": 180,
        "target": "Original H8 Hamiltonian/MPS, controls, four half-unit phases, and uncertainty budgets",
        "no_full_state_teacher": True, "independent_checker_process": True,
        "max_additional_observables": max_additional,
        "retention": "Every branch, including loose enclosures and resource failures; no discovery-improvement claim from this exploratory study",
        "sources": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "supplied_proposal": {"path": str(supplied), "sha256": hashlib.sha256(supplied.read_bytes()).hexdigest(),
                              "origin": read(supplied).get("recovery", {})} if supplied else None}
    write(output / "protocol.json", protocol)
    for p in files:
        archived = output / "source_snapshot" / p.name
        archived.parent.mkdir(exist_ok=True)
        archived.write_bytes(p.read_bytes())
    protocol_node = research.observe(action=Action.OPEN, state=State.OPEN,
        claim="Fixed independent teacher-free representation experiment", evidence={"protocol": protocol})

    def execute(action, request, label, timeout):
        request_file = output / (label + "_request.json")
        result_file = output / (label + ".json")
        write(request_file, {"inputs": str(copied), **request})
        with ledger.measure("frontier_" + action, label_id=label) as cost:
            process = run([sys.executable, "-B", "-S", "-m", MODULE, action,
                           str(request_file), str(result_file)], output / "processes" / label, timeout)
            cost.update(process)
            for name in ("stdout", "stderr"):
                cost[name + "_sha256"] = ledger.blob(Path(process[name]).read_bytes())
            if process["exit_code"] or not result_file.exists():
                cost["status"] = "timeout" if process["timeout"] else "failed"
                return {"status": State.COST.value if process["timeout"] else State.REPAIRABLE.value}
            result = read(result_file)
            cost["output_sha256"] = ledger.blob(result_file.read_bytes())
            cost["scientific_state"] = result["status"]
            return result

    baseline = execute("construct", {"source": BASE_SOURCE}, "baseline_observables", 120)
    baseline = execute("check", {"observables": str(output / "baseline_observables.json")}, "baseline_certificate", 180)
    baseline_node = research.observe(action=Action.PROVE, state=baseline["status"], reads=(protocol_node,),
        claim="Fixed four-observable baseline checked from original inputs", evidence=baseline)
    history, read_nodes = [], [protocol_node, baseline_node]
    for index in range(1 if supplied else 2):
        action = Action.OPEN if not history else choose_action(history[-1]["status"], [])
        prompt = PROMPT + "\nSelected research action: " + action.value + "\nPrior scoped outcomes:\n" + json.dumps([
            {k: value.get(k) for k in ("status", "D_interval", "error", "dimensions", "defect_norms")}
            for value in [baseline, *history]])
        if history and (output / f"proposal_{index-1}.json").exists():
            prompt += "\nPrevious program:\n" + read(output / f"proposal_{index-1}.json")["source"]
        if supplied:
            proposal = read(supplied)
            ledger.append("supplied_proposal_imported", source_sha256=ledger.blob(proposal["source"]),
                          provenance=protocol["supplied_proposal"], within_original_proposal_budget=False)
        else:
            proposal = propose(ledger, prompt, "representation_" + str(index), timeout=150)
        if proposal:
            write(output / f"proposal_{index}.json", proposal)
            label = "candidate_" + str(index)
            result = execute("construct", {"source": proposal["source"], "max_additional": max_additional}, label + "_observables", 120)
            if result["status"] == "proposal_constructed":
                request = {"observables": str(output / (label + "_observables.json"))}
                result = execute("check", request, label + "_certificate", 180)
                if result["status"] in (State.ACCEPTED.value, State.CONDITIONAL.value):
                    replayed = execute("replay", {**request, "certificate": str(output / (label + "_certificate.json"))},
                                       label + "_replay", 180)
                    if replayed != result:
                        result = {"status": State.REPAIRABLE.value, "error": "Fresh exact replay did not match", "replay": replayed}
        else:
            events = [row for row in ledger.events("operation_finished") if row["label"] == "model_proposal"]
            timed_out = bool(events and events[-1]["result"].get("timeout"))
            result = {"status": State.COST.value if timed_out else State.REPAIRABLE.value,
                      "error": "No executable native proposal within budget" if timed_out else "Native proposal unavailable"}
        node = research.observe(action=action, state=result["status"], reads=tuple(read_nodes),
            claim="Candidate observable representation under the unchanged physical task", evidence=result)
        if proposal:
            method = Method("observable_representation_" + str(index), ledger.blob(proposal["source"]),
                {"payload": "original H8 Gaussian-rational operators; no teacher"},
                {"observables": "additional exact Hermitian particle-neutral polynomials"},
                (proposal["assumptions"], "Every new output requires independent defect checking"),
                ("control-observable-representation",), (node,), {"budget": protocol}, tuple(read_nodes),
                (proposal["failure_cases"], result.get("remaining_obligation", result.get("error", ""))),
                "propose(payload); then independent fixed-task construct/replay", 
                "enabling" if result["status"] == State.ACCEPTED.value else "high_risk")
            research.retain(method)
        read_nodes.append(node)
        history.append(result)
        print(json.dumps({"candidate": index, "status": result["status"],
                          "interval": result.get("D_interval"), "dimensions": result.get("dimensions")}), flush=True)
    result = {"baseline": {k: baseline.get(k) for k in ("status", "D_interval", "dimensions")},
        "proposals": [{k: r.get(k) for k in ("status", "D_interval", "dimensions", "error", "remaining_obligation")} for r in history],
        "target_proved": any(r.get("target_proved", False) for r in history),
        "supplied_proposal": protocol["supplied_proposal"],
        "discovery_improvement_proved": False, "costs": ledger.costs(), "audit": ledger.audit()}
    write(output / "result.json", result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    action = sys.argv[1]
    if action == "campaign":
        campaign(Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve(),
                 Path(sys.argv[4]).resolve() if len(sys.argv) > 4 else None,
                 int(sys.argv[5]) if len(sys.argv) > 5 else 14)
    else:
        worker(action, *(Path(x).resolve() for x in sys.argv[2:]))
