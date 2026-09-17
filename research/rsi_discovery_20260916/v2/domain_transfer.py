"""Later-process reactivation on all actual SOS operator channels.

This is an additional conformance gate, not another independent campaign and
not a replacement for the frozen four-arm primary test.
"""
import json
from pathlib import Path
import sys
import time

from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.process import run
from research.rsi_discovery_20260916.v2.algebra import oracle, gram_reference, validate_map
from research.rsi_discovery_20260916.v2.library import MethodLibrary


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def worker(campaign, arm, payload):
    result = json.loads((campaign / "result.json").read_text())["arms"][arm]
    if not result.get("method_id"):
        print(json.dumps({"status": "unobserved", "reason": "No newly discovered method was retained"})); return
    started = time.perf_counter()
    library = MethodLibrary(Ledger(campaign))
    method = library.activate(result["method_id"], {"oracle": oracle})
    activation = time.perf_counter() - started
    records = []
    for group in json.loads(payload.read_text()):
        words = [tuple(tuple(letter) for letter in word) for word in group["words"]]
        start = time.perf_counter()
        try:
            proposed = method({"words": list(words), "cache": {}})
            validate_map(proposed, words)
            expected = gram_reference(words)
            if proposed != expected:
                raise ValueError("A retained map disagrees with the original exact coefficient semantics")
        except Exception as error:
            print(json.dumps({"status": "implementation_failure", "group": group["name"],
                "error": str(error), "completed_groups": len(records)})); return
        records.append({"name": group["name"], "dictionary_words": len(words),
                        "word_charge": sum(2 * c - 1 for c, _ in words[0]),
                        "max_word_degree": max(map(len, words)), "seconds": time.perf_counter() - start})
    if any(k in sys.modules for k in ("numpy", "scipy", "cvxpy", "pyscf", "quimb")):
        raise AssertionError("Numerical imports in the later-process exact gate")
    print(json.dumps({"status": "verified_encoded_claim", "activation_seconds": activation,
        "method_id": result["method_id"], "source_sha256": result["frozen_source"], "groups": records,
        "new_model_calls": 0, "method_changed_after_freezing": False, "new_Hamiltonian_quality_claim": False}))


def main(study, output):
    if not (study / "result.json").exists():
        raise ValueError("Finish the frozen campaign before later-process conformance")
    ledger = Ledger(output)
    protocol = {"type": "additional_semantic_conformance", "modes": [8, 12],
        "clusters": "all consecutive width-three spatial clusters", "collective_pairs": True,
        "channels": "every emitted linear, pair, particle-hole, mixed and triple charge/spin group",
        "scientific_test": "exact raw symmetric map equality; unchanged source; original later-process dependency activation",
        "not_an_independent_RSI_campaign": True, "no_model_calls": True}
    ledger.append("conformance_protocol", protocol=protocol)
    with ledger.measure("generate_actual_operator_dictionaries"):
        from research.interacting_scaling_20260915.dictionary import clustered_frame
        payload = []
        for modes in protocol["modes"]:
            groups, _ = clustered_frame(modes, [tuple(range(i, i + 3)) for i in range(modes // 2 - 2)])
            payload.extend({"name": str(modes) + ":" + g["name"], "words": g["words"]} for g in groups if g["words"])
        write(output / "payload.json", payload)
        ledger.append("dictionary_payload_frozen", sha256=ledger.blob((output / "payload.json").read_bytes()))
    results = []
    for campaign in sorted(study.glob("campaign_*")):
        for arm in "ABCD":
            label = campaign.name + "_" + arm
            with ledger.measure("later_process_exact_transfer", campaign=campaign.name, arm=arm) as cost:
                process = run([sys.executable, "-B", "-S", "-m", "research.rsi_discovery_20260916.v2.domain_transfer",
                    "worker", str(campaign), arm, str(output / "payload.json")], output / label, 120)
                cost.update(process)
                value = json.loads(Path(process["stdout"]).read_text()) if not process["exit_code"] else {
                    "status": "resource_obstruction" if process["timeout"] else "implementation_failure"}
                cost["scientific_state"] = value["status"]
                cost["stdout_sha256"] = ledger.blob(Path(process["stdout"]).read_bytes())
            write(output / (label + ".json"), value)
            results.append({"campaign": campaign.name, "arm": arm, "status": value["status"],
                            "groups": len(value.get("groups", [])), "activation_seconds": value.get("activation_seconds")})
            print(json.dumps(results[-1]), flush=True)
    write(output / "result.json", {"protocol": protocol, "results": results,
                                   "costs": ledger.costs(), "audit": ledger.audit()})


if __name__ == "__main__":
    if sys.argv[1] == "worker":
        worker(Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]))
    else:
        main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
