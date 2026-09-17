"""Structural correlation proposals consumed by the actual non-enumerating solver."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

from experiments.marginal_symbolic import decode
from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.process import ROOT, run
from research.rsi_discovery_20260916.v2.campaign import write, read
from research.rsi_discovery_20260916.v2.language import compile_program
from research.rsi_discovery_20260916.v2.model import propose
from research.rsi_discovery_20260916.v2.records import Research, Method, State, Action, environment


def validate_design(value, spatial):
    if type(value) is not dict or set(value) != {"clusters", "collective_clusters"}:
        raise ValueError("A representation may choose only local and collective supports")
    result, seen, count = {}, set(), 0
    for name in ("clusters", "collective_clusters"):
        if type(value[name]) is not list:
            raise ValueError("Support collections must be lists")
        result[name] = []
        for support in value[name]:
            if type(support) not in (list, tuple) or not 2 <= len(support) <= min(3, spatial):
                raise ValueError("Use direct two- or three-orbital supports")
            if any(type(i) is not int or not 0 <= i < spatial for i in support) or len(set(support)) != len(support):
                raise ValueError("Invalid support in the original orbital domain")
            key = tuple(sorted(support))
            if key in seen:
                raise ValueError("Duplicate support wastes the declared search budget")
            seen.add(key); count += 1
            result[name].append(list(key))
    if not 1 <= count <= max(1, spatial - 2):
        raise ValueError("Support count exceeds the matched width-three baseline budget")
    return result


def select(data, source):
    if data["modes"] % 2 or data["modes"] < 4:
        raise ValueError("Paired spin orbitals are required")
    spatial = data["modes"] // 2
    # No Gram map, sector states or teacher moments exist when this program runs.
    payload = {"spatial": spatial, "modes": data["modes"], "particles": data["particles"],
        "hamiltonian": decode(data["hamiltonian"], data["modes"], 4),
        "max_supports": max(1, spatial - 2), "max_support_size": min(3, spatial), "cache": {}}
    return validate_design(compile_program(source, {"rational": F})(payload), spatial)


PROMPT = """Discover an executable Spectra procedure for selecting useful collective correlations
directly from the original molecular Hamiltonian, before constructing any global cubic Gram map.
Return pure Python propose(payload). Input: spatial (number of spatial orbitals), modes=2*spatial,
particles, hamiltonian (canonical CAR word -> exact Fraction), max_supports=max(1,spatial-2),
max_support_size=min(3,spatial), and an empty cache. Spatial orbital p has spin modes 2p and 2p+1.
An exact scalar coefficient factory rational(n,d) is available. No state, occupation list,
ground-state vector, full coefficient map, or trajectory is supplied.
Output exactly {'clusters': [...], 'collective_clusters': [...]}. Each support is a unique
list of two or three distinct spatial indices. Across both lists, use 1..max_supports supports.
The existing global quadratic and collective pair-supported channels remain included.
Each 'clusters' entry creates its own cubic Gram channels. All 'collective_clusters' supports
join the shared collective pair Gram channels, retaining cross terms between distant supports.
The current number-ideal multiplier dictionary uses the local 'clusters' list; moving a support
into 'collective_clusters' expands shared PSD coupling but does not add that support's local
number-ideal multipliers. Account for that actual constructor tradeoff.
You may use nonlocal/overlapping supports, weighted integral graphs, a different partition, or
another structural selection algorithm. Do not merely change a numeric solver parameter.
The frozen comparator uses consecutive width-three local clusters and no collective clusters.
Your procedure must work at new orbital counts. A four-H nonuniform STO-3G molecule is the
development case. A larger molecular input is withheld until the source is frozen.
The downstream workflow constructs fresh MPS and exact SOS bounds without full-sector enumeration.
Both spin sectors receive the same chosen supports; the independent checker accepts an energy
interval against the original H and full fixed-N sector, with chemical target 1.6 mHa.
A negative moment or fewer supports is not success. More collective coupling may increase map,
Gram and solve cost. The actual two-million-entry Gram and 180000-row caps still apply.
Supply one concise general structural method. No external tools, files or model calls.
"""


def campaign(output, map_proposal, primitive_proposal):
    output.mkdir(parents=True, exist_ok=False)
    ledger = Ledger(output)
    files = [Path(__file__), Path(__file__).with_name("workflow.py"), Path(__file__).with_name("language.py")]
    baseline_root = ROOT / "results/rsi_discovery_v2_20260916"
    def original_specification(name):
        target = read(baseline_root / name / "result.json")["target"]
        return {**target["specification"], "particles": target["N"]}
    development = original_specification("nonenumerating_reference001")
    transfer = original_specification("physical_transfer_h6001")
    protocol = {"native_model": "gpt-6-astra", "effort": "high", "proposal_calls": 1, "model_deadline_seconds": 300,
        "development": development, "transfer_hash": digest(transfer), "no_repair_from_transfer": True,
        "budgets": {"development": {"bond": 16, "solve_seconds_per_sector": 50},
                    "transfer": {"bond": 24, "solve_seconds_per_sector": 60}},
        "comparison_scope": "One exploratory structural algorithm against previously measured frozen windows; not an independent RSI significance test",
        "map_source_sha256": digest(read(map_proposal)["source"]),
        "primitive_source_sha256": digest(read(primitive_proposal)["source"]),
        "sources": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    write(output / "protocol.json", protocol)
    for path in files:
        archived = output / "source_snapshot" / path.name
        archived.parent.mkdir(exist_ok=True); archived.write_bytes(path.read_bytes())
    research = Research(ledger, {"claim": "Useful directly selected collective structure at original chemical tolerance"},
                        environment({"development": development, "transfer_hash": digest(transfer)}, files,
                                    [protocol["map_source_sha256"], protocol["primitive_source_sha256"]]))
    origin = research.observe(action=Action.OPEN, state=State.OPEN,
        claim="Same original target; only representation supports may change", evidence={"protocol": protocol})
    proposal = propose(ledger, PROMPT, "select-correlations", timeout=300)
    if not proposal:
        record = ledger.events("operation_finished")[-1]["result"]
        write(output / "result.json", {"status": State.COST.value if record.get("timeout") else State.REPAIRABLE.value,
              "reason": "No executable native proposal",
              "discovery_improvement_proved": False, "costs": ledger.costs(), "audit": ledger.audit()})
        return
    proposal_path = output / "proposal.json"; write(proposal_path, proposal)
    ledger.append("source_frozen_before_transfer", source_sha256=ledger.blob(proposal["source"]))
    outcomes = []
    for label, specification, bond, seconds in (("development", development, 16, 50), ("transfer", transfer, 24, 60)):
        spec_path = output / (label + "_specification.json");write(spec_path, specification)
        case = output / label
        # The workflow owns all nested stage costs. Do not wrap it in another
        # metered operation and double-charge the same model/solver wall time.
        process = run([sys.executable, "-B", "-m", "research.rsi_discovery_20260916.v2.workflow", "campaign", str(case),
            str(map_proposal), str(primitive_proposal), "--specification", str(spec_path), "--bond", str(bond),
            "--seconds", str(seconds), "--selector", str(proposal_path)], output / "processes" / label, 1500)
        ledger.append("workflow_dispatch_completed", label=label, process=process,
                      accounting="Nested workflow stage operations are charged in their own ledger")
        value = read(case / "result.json") if (case / "result.json").exists() else {
            "status": State.COST.value if process["timeout"] else State.REPAIRABLE.value, "acceptance": {}}
        current = environment({"case": label, "specification": specification, "bond": bond, "solve_seconds": seconds},
                              files, [digest(proposal["source"]), protocol["map_source_sha256"], protocol["primitive_source_sha256"]])
        research.env = current["fingerprint"]
        ledger.append("live_environment_transition", environment=current)
        node = research.observe(action=Action.TRANSFER, state=value["status"], reads=(origin,),
            claim="Original full-N molecular interval using directly selected supports", evidence=value)
        outcomes.append({"case": label, "status": value["status"], "acceptance": value["acceptance"],
                         "costs": value.get("costs"), "node": node})
        print(json.dumps({"case": label, "status": value["status"], "width_mHa": value["acceptance"].get("width_mHa")}), flush=True)
    accepted = [r["node"] for r in outcomes if r["status"] == State.ACCEPTED.value]
    method = Method("direct_correlation_supports", ledger.blob(proposal["source"]),
        {"payload": "original Hamiltonian and orbital domain, no global map or teacher"},
        {"design": "local and coupled collective supports"},
        (proposal["assumptions"], "Every new physical input requires independent full-N acceptance; scope is the recorded inputs only"),
        ("molecular-correlation-selection",), tuple(accepted or [r["node"] for r in outcomes]),
        {"budget": protocol["budgets"], "outcomes": outcomes}, (origin,),
        (proposal["failure_cases"], "A legal support design does not imply chemical accuracy or lower complete cost"),
        "select(original_fixture, source); then exact full-N workflow acceptance", "enabling" if accepted else "high_risk")
    retained = research.retain(method)
    write(output / "result.json", {"outcomes": outcomes, "retained_method": retained,
        "discovery_improvement_proved": False, "costs": ledger.costs(), "audit": ledger.audit(),
        "scope": protocol["comparison_scope"]})


if __name__ == "__main__":
    if sys.argv[1] == "worker":
        fixture, source, output = map(Path, sys.argv[2:])
        write(output, select(read(fixture), read(source)["source"]))
    else:
        campaign(*(Path(x).resolve() for x in sys.argv[2:]))
