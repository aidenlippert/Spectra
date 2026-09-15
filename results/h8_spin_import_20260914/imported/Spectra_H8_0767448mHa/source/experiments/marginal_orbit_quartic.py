"""Symmetry-average quartic CAR certificates for the matched hopping model.

The averaging is exact at the operator level: each Gram matrix is projected
onto the commutant of the S_m x C_2 flavor symmetry.  It therefore cannot
decrease the SOS bound for an invariant Hamiltonian, and the resulting
rounded certificate is checked by the existing exact verifier.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json
import numpy as np

from experiments.marginal_coefficient import (dictionaries, hopping_model,
    solve_coefficients, symmetric_upper, export)
from experiments.marginal_symbolic import mono, transform, add, scale


def symmetry_group(modes):
    m = modes // 2
    if modes % 2:
        raise ValueError("matched model requires even mode count")
    return [tuple(p[i % m] + m * ((i // m) ^ swap)
                  for i in range(modes))
            for p in permutations(range(m)) for swap in (0, 1)]


def word_action(word, mapping):
    image = transform(mono(word), mapping)
    if len(image) != 1:
        raise ValueError("symmetry sent a dictionary word to zero or a sum")
    return next(iter(image.items()))


def average_solution(solution, blocks, modes):
    """Project numerical Gram matrices to the exact permutation commutant."""
    group = symmetry_group(modes)
    out = []
    for block, q in zip(blocks, solution["grams"]):
        # Raw dictionary words need not be in the same annihilator ordering
        # after a permutation.  Index their canonical CAR polynomials.
        index = {}
        for i, w in enumerate(block["words"]):
            key = next(iter(transform(mono(w), tuple(range(modes)))))
            index[key] = i
        acc = np.zeros_like(q)
        for mapping in group:
            p = np.zeros_like(q)
            for i, word in enumerate(block["words"]):
                image, sign = word_action(word, mapping)
                key = next(iter(transform(mono(image), tuple(range(modes)))))
                j = index.get(key)
                if j is None:
                    raise ValueError(f"dictionary is not symmetry closed: {image}")
                p[j, i] = sign
            acc += p.T @ q @ p
        out.append((acc / len(group) + (acc / len(group)).T) / 2)
    result = dict(solution)
    result["grams"] = out
    # The number multiplier must be projected by the same group action.
    if solution["basis"]:
        x = add(*(scale(q, float(a)) for q, a in zip(solution["basis"], solution["x"])))
        xp = {}
        for mapping in group:
            y = transform(x, mapping)
            for w, a in y.items(): xp[w] = xp.get(w, 0) + a / len(group)
        words = sorted(set().union(*(set(q) for q in solution["basis"])), key=lambda w: (len(w), w))
        A = np.array([[float(q.get(w, 0)) for q in solution["basis"]] for w in words])
        result["x"] = np.linalg.lstsq(A, np.array([float(xp.get(w, 0)) for w in words]), rcond=None)[0]
    return result


def run(modes=6, t=F(1, 5), family="mixed"):
    h = hopping_model(modes, t)
    blocks = dictionaries(modes, family)
    raw = solve_coefficients(h, modes, modes // 2, blocks)
    projected = average_solution(raw, blocks, modes)
    certificate, receipt = export(h, modes, modes // 2, blocks, projected,
                                  operator_degree=4)
    certificate["symmetry"] = {"group": "S_m x C2", "order": len(symmetry_group(modes)),
                                "m": modes // 2}
    upper = symmetric_upper(modes, t)
    certificate["variational_upper"] = upper
    receipt.update({"modes": modes, "t": str(t), "family": family,
                    "symmetry_order": len(symmetry_group(modes)),
                    "upper": upper["upper"],
                    "width_float": float(F(upper["upper"]) - F(receipt["lower"]))})
    return certificate, receipt


if __name__ == "__main__":
    c, r = run()
    out = Path(__file__).resolve().parents[1] / "results/marginal_orbit_quartic"
    out.mkdir(exist_ok=True)
    (out / "m6_t1_5_mixed.json").write_text(json.dumps(c) + "\n")
    (out / "m6_t1_5_mixed_receipt.json").write_text(json.dumps(r, indent=2) + "\n")
    print(json.dumps(r), flush=True)
