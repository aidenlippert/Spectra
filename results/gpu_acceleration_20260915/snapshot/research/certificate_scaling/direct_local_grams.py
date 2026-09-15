"""Bounded SOS proposals whose Gram dictionaries come only from H locality.

The Hamiltonian support is scanned for mode hyperedges.  Each hyperedge gets
small same-charge annihilation/creation blocks built from subsets of that
edge; no source certificate factors or supports are read.  Numerical CVXPY
proposals are accepted only after exact CAR replay.
"""
from __future__ import annotations

import argparse, json, sys, time
from itertools import combinations
from pathlib import Path
from fractions import Fraction as F

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_coefficient import solve_coefficients, export
from experiments.marginal_symbolic import decode


def dagger(word):
    return tuple((1 - c, i) for c, i in reversed(word))


def locality_blocks(h, modes: int, max_words: int = 4):
    edges = set()
    for word in h:
        edge = tuple(sorted({i for _, i in word}))
        if edge:
            edges.add(edge)
    groups = {}
    for edge in sorted(edges, key=lambda x: (len(x), x)):
        for degree in range(0, min(3, len(edge)) + 1):
            for subset in combinations(edge, degree):
                w = tuple((0, i) for i in subset)
                groups.setdefault((edge, -degree), []).append(w)
                groups.setdefault((edge, degree), []).append(dagger(w))
    blocks = []
    for (edge, charge), words in groups.items():
        words = list(dict.fromkeys(words))
        if not words or charge == 0:
            continue
        words = words[:max_words]
        blocks.append({"name": f"local-edge-{edge}-q{charge}", "words": words})
    return blocks


def run(input_path: Path, output: Path, receipt_path: Path, max_words: int, denominator: int):
    data = json.loads(input_path.read_text())
    modes = data["modes"]
    particles = data.get("particles", modes // 2)
    h = decode(data["hamiltonian"], modes, 4)
    blocks = locality_blocks(h, modes, max_words)
    started = time.monotonic()
    solution = solve_coefficients(h, modes, particles, blocks, True)
    certificate, receipt = export(h, modes, particles, blocks, solution, denominator)
    receipt.update({"input": str(input_path), "max_words": max_words,
                    "local_block_count": len(blocks),
                    "local_block_max_dimension": max((len(b["words"]) for b in blocks), default=0),
                    "local_gram_scalar_variables": sum(len(b["words"])**2 for b in blocks),
                    "wall_seconds": time.monotonic() - started,
                    "source_upper_candidates": {k: data[k] for k in ("independent_upper", "car_rational_ground", "total_fci") if k in data}})
    output.parent.mkdir(parents=True, exist_ok=True); receipt_path.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(certificate, indent=2) + "\n")
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--receipt", type=Path, required=True)
    p.add_argument("--max-words", type=int, default=4)
    p.add_argument("--denominator", type=int, default=10**7)
    args = p.parse_args()
    print(json.dumps(run(args.input, args.output, args.receipt, args.max_words, args.denominator), indent=2))


if __name__ == "__main__": main()
