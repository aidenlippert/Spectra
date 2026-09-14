"""Discovery and exact replay driver for the joint half/charged overlap bound."""
import hashlib, json, pathlib, sys, time
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_charged_projectors import (  # noqa: E402
    charged_projector_bound, joint_overlap_grams, joint_projector_bound,
)

OUT = ROOT / "results/marginal_graded_hubbard8/joint_overlap_exact"
CHARGED = ROOT / "results/marginal_graded_hubbard8/charged_projector/source.json"
HALF = ROOT / "results/marginal_graded_hubbard8/six_site_projector/refined_certificate.json"

def load_sources():
    charged = json.loads(CHARGED.read_text())["physical_state"]
    half = json.loads(HALF.read_text())["vector"]
    return {int(k): int(v) for k, v in half.items()}, {int(k): int(v) for k, v in charged.items()}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    half, charged = load_sources()
    receipt = {"driver": sha(pathlib.Path(__file__)), "sources": {str(CHARGED): sha(CHARGED), str(HALF): sha(HALF)}, "started_unix": time.time()}
    t = time.perf_counter()
    charged_result = charged_projector_bound(charged, 4, Fraction(264879339, 125000000))
    receipt["charged_replay"] = charged_result
    receipt["charged_seconds"] = time.perf_counter() - t
    groups = joint_overlap_grams(half, charged, 4)
    receipt["joint_counts"] = {"columns": sum(len(x["gram"]) for x in groups.values()), "max_sector": max(len(x["gram"]) for x in groups.values()), "sectors": len(groups)}
    receipt["joint_proposals"] = []
    for ratio in (Fraction(1, 2), Fraction(1), Fraction(2)):
        # Conservative outward ceilings are supplied by discovery; these are exact replay candidates.
        ceiling = {Fraction(1,2): Fraction(595736407,250000000), Fraction(1): Fraction(189969857,62500000), Fraction(2): Fraction(4925489283,1000000000)}[ratio]
        t = time.perf_counter()
        result = joint_projector_bound(half, charged, 4, ratio, ceiling)
        result["seconds"] = time.perf_counter() - t
        receipt["joint_proposals"].append(result)
    receipt["finished_unix"] = time.time()
    (OUT / "joint_overlap_replay.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"columns": receipt["joint_counts"], "charged": receipt["charged_replay"]["average_fidelity_ceiling"], "proposals": [(x["charged_ratio"], x["average_ceiling"]) for x in receipt["joint_proposals"]]}, indent=2))

if __name__ == "__main__":
    main()
