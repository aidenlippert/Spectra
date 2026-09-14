"""Bounded sparse/chordal certificate probe on the canonical H4 fixture.

This is deliberately an adapter around the repository's exact CAR verifier:
the experiment measures the support graph induced by the adaptive quadratic
blocks, then replays the certificate and exercises a malformed-factor refusal.
It does not claim that the sparse family is globally optimal.
"""
from __future__ import annotations

import copy, json, time, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from experiments.marginal_symbolic import verify

SOURCE = ROOT / "results/certificate_scaling/active_space_ladder/h4/fixture.json"


def probe() -> dict:
    started = time.monotonic()
    fixture = json.loads(SOURCE.read_text())
    from research.certificate_scaling.direct_sparse_discovery import discover
    from experiments.marginal_symbolic import decode
    h = decode(fixture["hamiltonian"], fixture["modes"], 4)
    cert, exact = discover(h, fixture["modes"], fixture["particles"], budget=64, time_limit=30)
    blocks = cert["blocks"]
    words = [tuple(tuple(x) for x in w) for b in blocks for w in b["words"]]
    vertices = {w for w in words}
    edges = set()
    for b in blocks:
        ws = [tuple(tuple(x) for x in w) for w in b["words"]]
        for i, left in enumerate(ws):
            for right in ws[i + 1:]:
                edges.add(tuple(sorted((left, right), key=repr)))
    # Chordal-style local blocks: treewidth is bounded by the largest block
    # support minus one; this is a structural cost upper bound, not a theorem.
    malformed = copy.deepcopy(cert)
    malformed["blocks"][0]["factor"][0][0] = "corrupt"
    try:
        verify(malformed)
    except (TypeError, ValueError, KeyError, ZeroDivisionError) as exc:
        refusal = {"accepted": False, "error": type(exc).__name__}
    else:
        refusal = {"accepted": True, "error": None}
    return {
        "hypothesis": "Adaptive local quadratic blocks retain a sparse interaction graph while exact replay preserves the lower bound.",
        "falsifier": "A malformed factor must be rejected; a future H6/H10 run falsifies scaling if block support or full pricing grows superlinearly.",
        "source": str(SOURCE), "method": "Hamiltonian-only direct H-ranked sparse quadratic LP (64-atom family)",
        "baseline_provenance": "The paired 512 run is a larger512-atomLP, not the full quadratic PSD/SOS cone.",
        "blocks": len(blocks), "block_word_count": len(words), "unique_words": len(vertices),
        "interaction_edges": len(edges),
        "max_block_order": max(len(b["words"]) for b in blocks),
        "factor_nonzeros": sum(sum(1 for x in row if x) for b in blocks for row in b["factor"]),
        "exact_replay": {"lower_float": exact["lower_float"], "residual_l1_float": exact["residual_l1_float"]},
        "corrupt_replay": refusal, "wall_seconds": time.monotonic() - started,
        "scope": "Frozen rational H4 Hamiltonian; discovery cost and omitted-family optimality are not certified.",
    }


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2))
