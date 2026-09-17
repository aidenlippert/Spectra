"""Connect retained programs to the real local/collective SOS + charge-MPS workflow.

All accepting paths run in separate standard-library-only processes. The proposed
maps change construction, never the Hamiltonian, sector, tolerance or checker.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import shutil
import sys
import time

from research.rsi_discovery_20260916.ledger import Ledger, digest
from research.rsi_discovery_20260916.process import run


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def prepare(case, method, primitive):
    import numpy as np
    from scipy import sparse
    from research.interacting_scaling_20260915 import prepare as current
    from research.rsi_discovery_20260916.v2.algebra import oracle, validate_map
    from research.rsi_discovery_20260916.v2.language import compile_program
    capabilities = {"oracle": oracle}
    if primitive:
        capabilities["car_kernel"] = compile_program(primitive, capabilities)
    candidate = compile_program(method, capabilities)
    receipts = []
    original = current.gram_map

    def proposed_map(words, lookup):
        started = time.perf_counter()
        # Each activation pays a fresh cache cost, as in the discovery/transfer gate.
        value = candidate({"words": words, "cache": {}})
        validate_map(value, words)
        row, col, values = [], [], []
        size = len(words)
        for word, entries in value.items():
            if word not in lookup:
                continue
            for (i, j), c in entries.items():
                row.append(lookup[word]); col.append(i * size + j)
                values.append(c if i == j else c / 2)
                if i != j:
                    row.append(lookup[word]); col.append(j * size + i); values.append(c / 2)
        matrix = sparse.csr_matrix((values, (row, col)), shape=(len(lookup), size * size))
        construction = time.perf_counter() - started
        # Exact half-integer arrays: this comparison introduces no rounding tolerance.
        check_start = time.perf_counter()
        reference = original(words, lookup)
        order = (np.arange(size)[None, :] * size + np.arange(size)[:, None]).ravel()
        reference = (reference + reference[:, order]) * .5
        difference = matrix - reference
        difference.eliminate_zeros()
        if difference.nnz:
            raise ValueError("Proposed map changed an actual solver coefficient")
        receipts.append({"words": size, "construction_seconds": construction,
                         "independent_check_seconds": time.perf_counter() - check_start,
                         "nonzeros": matrix.nnz, "exact_symmetric_map_equal": True})
        return matrix

    current.gram_map = proposed_map
    current.build(case)
    write(case / "prepared" / "method_integration.json", {"source_sha256": digest(method),
        "primitive_sha256": digest(primitive), "blocks": receipts,
        "original_semantics_preserved": True, "energy_coordinate_not_in_ideal": True})


def campaign(output, method_path, primitive_path, specification=None, bond=16, seconds=50):
    output = output.resolve()
    ledger = Ledger(output)
    method = json.loads(method_path.read_text())["source"]
    primitive = json.loads(primitive_path.read_text())["source"]
    specification = specification or {"name": "h4_nonuniform_v2", "basis": "sto-3g", "particles": 4,
        "geometry": [["H", [0., 0., z]] for z in (0., 1.05, 2.31, 3.62)]}
    input_dir = output / "fresh_input"
    input_dir.mkdir()
    write(input_dir / "specification.json", specification)
    target = {"specification": json.loads((input_dir / "specification.json").read_text()),
              "N": specification["particles"], "scope": "entire fixed-N sector", "max_width_Ha": "1/625",
              "constructor": "declared collective pair supports plus width-three local clusters",
              "full_occupation_enumeration": False, "cached_teacher": False}
    ledger.append("target_frozen", target=target, sha256=digest(target))

    def command(label, arguments, seconds=120, stdlib=False):
        with ledger.measure(label) as cost:
            result = run([sys.executable, "-B", *(["-S"] if stdlib else []), "-m", *arguments],
                         output / "processes" / label, seconds)
            cost.update(result)
            for stream in ("stdout", "stderr"):
                cost[stream + "_sha256"] = ledger.blob(Path(result[stream]).read_bytes())
            if result["exit_code"]:
                cost["status"] = "timeout" if result["timeout"] else "failed"
                raise RuntimeError(label + " failed: " + result["stderr"])
        return result

    try:
        command("fresh_integrals", ["research.transfer_solver_20260915.generate", str(input_dir)], 90)
        fixture = json.loads((input_dir / "fixture.json").read_text())
        if fixture["particles"] != specification["particles"]:
            raise ValueError("Generated particle number differs from frozen target")
        spatial = fixture["modes"] // 2
        width = min(3, spatial)
        clusters = [list(range(i, i + width)) for i in range(spatial - width + 1)]
        command("fresh_MPS", ["research.interacting_scaling_20260915.state", str(input_dir), "--bond", str(bond), "--sweeps", "5"], 150)
        command("exact_MPS_upper", ["research.correlated_pair_20260913.mps_exact", str(input_dir / "fixture.json"),
                str(input_dir / "mps/state.json"), "--output", str(input_dir / "upper.json")], 90, True)
        for mag in (0, 1):
            case = output / ("singlet" if mag == 0 else "nonsinglet")
            case.mkdir()
            for name in ("fixture.json", "upper.json"):
                shutil.copyfile(input_dir / name, case / name)
            (case / "mps").mkdir()
            shutil.copyfile(input_dir / "mps/state.json", case / "mps/state.json")
            write(case / "design.json", {"clusters": clusters, "collective_pairs": True,
                "complete": False, "magnetization": mag, "max_Gram_entries": 2000000,
                "max_coefficient_rows": 180000, "global_particle_number_only": True})
            command("prepare_" + str(mag), ["research.rsi_discovery_20260916.v2.workflow", "prepare", str(case), str(method_path.resolve()), str(primitive_path.resolve())], 150)
            command("solve_" + str(mag), ["research.interacting_scaling_20260915.solve", str(case), "discovered",
                    "--seconds", str(seconds), "--mu", "2"], seconds + 90)
        command("independent_full_N_acceptance", ["research.collective_completion_20260914.spin_screen",
            str(input_dir / "fixture.json"), str(output / "singlet/discovered/export/certificate.json"),
            str(output / "nonsinglet/discovered/export/certificate.json"), str(output / "acceptance.json"),
            "--upper", str(input_dir / "upper.json")], 150, True)
        acceptance = json.loads((output / "acceptance.json").read_text())
        status = "verified_encoded_claim" if acceptance["target_1p6mHa_met"] else "conditional_proof_obligation"
        if not acceptance["target_1p6mHa_met"]:
            acceptance["remaining_obligation"] = "Interval exceeds target after the fixed optimization budget; representational impossibility is not proved."
    except Exception as error:
        acceptance, status = {"error": str(error)}, "implementation_failure"
    result = {"status": status, "target": target, "acceptance": acceptance,
              "costs": ledger.costs(), "audit": ledger.audit()}
    write(output / "result.json", result)
    print(json.dumps(result))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=("campaign", "prepare"))
    p.add_argument("output", type=Path)
    p.add_argument("method", type=Path)
    p.add_argument("primitive", type=Path)
    p.add_argument("--specification", type=Path)
    p.add_argument("--bond", type=int, default=16)
    p.add_argument("--seconds", type=int, default=50)
    a = p.parse_args()
    if a.action == "prepare":
        prepare(a.output, json.loads(a.method.read_text())["source"], json.loads(a.primitive.read_text())["source"])
    else:
        campaign(a.output, a.method, a.primitive,
                 json.loads(a.specification.read_text()) if a.specification else None, a.bond, a.seconds)


if __name__ == "__main__":
    main()
