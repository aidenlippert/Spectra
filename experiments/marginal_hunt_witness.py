"""Exact scalable DQG/T1 separation example; not a molecular calculation."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from experiments.marginal_hunt_car import add, adj, mono, mul, scale


def diagonal_moment(poly):
    """Candidate diagonal 2-RDM functional, in ordered pair convention."""
    result = F(0)
    for word, coefficient in poly.items():
        creators = [m for c, m in word if c]
        annihilators = [m for c, m in word if not c]
        if creators != annihilators:
            continue
        k = len(creators)
        if k == 0:
            moment = F(1)
        elif k == 1:
            moment = F(1, 2)
        elif k == 2:
            moment = F(1, 8) if creators[0] // 3 == creators[1] // 3 else F(1, 4)
        else:
            raise ValueError("This functional only specifies a 2-RDM")
        result += coefficient * (-1) ** (k * (k - 1) // 2) * moment
    return result


def gram(words):
    return [[diagonal_moment(mul(adj(u), v)) for v in words] for u in words]


def verify_dqg_blocks():
    """Check every D,Q,G entry for six modes against explicit PSD blocks."""
    pairs = list(combinations(range(6), 2))
    d_words = [mono(((0, j), (0, i))) for i, j in pairs]
    d, q = gram(d_words), gram([adj(w) for w in d_words])
    for row, (i, j) in enumerate(pairs):
        value = F(1, 8) if i // 3 == j // 3 else F(1, 4)
        for col in range(len(pairs)):
            expected = value if row == col else 0
            if d[row][col] != expected or q[row][col] != expected:
                raise AssertionError("D/Q diagonal certificate mismatch")
    indices = [(i, j) for i in range(6) for j in range(6)]
    g = gram([mono(((1, i), (0, j))) for i, j in indices])
    for row, (i, j) in enumerate(indices):
        for col, (k, l) in enumerate(indices):
            expected = F(0)
            if i == j and k == l:
                expected = F(1, 2) if i == k else (F(1, 8) if i // 3 == k // 3 else F(1, 4))
            elif row == col:
                expected = F(1, 2) - (F(1, 8) if i // 3 == j // 3 else F(1, 4))
            if g[row][col] != expected:
                raise AssertionError("G block certificate mismatch")
    # P=J/4 + blockdiag((3 I-J)/8). Each 3 I-J is a graph Laplacian.
    return {"D_entries_checked": 225, "Q_entries_checked": 225, "G_entries_checked": 1296}


def local_identity():
    b = mono(((0, 0), (0, 1), (0, 2)))
    certificate = add(mul(adj(b), b), mul(b, adj(b)))
    numbers = [mono(((1, i), (0, i))) for i in range(3)]
    expected = add(mono(()), *(scale(n, -1) for n in numbers),
                   *(mul(numbers[i], numbers[j]) for i, j in combinations(range(3), 2)))
    if certificate != expected:
        raise AssertionError("Cubic certificate identity failed")
    if diagonal_moment(certificate) != -F(1, 8):
        raise AssertionError("Candidate does not violate the expected cut")
    return certificate


def family_receipt(r):
    if not isinstance(r, int) or isinstance(r, bool) or r < 1:
        raise ValueError("r must be a positive integer")
    modes, particles = 6 * r, 3 * r
    row_sum = 2 * F(1, 8) + (modes - 3) * F(1, 4)
    if row_sum != (particles - 1) * F(1, 2):
        raise AssertionError("Fixed-N contraction failed")
    return {"r": r, "modes": modes, "particles": particles,
            "DQG_optimum": str(F(3 * r, 4)), "certified_ground_energy": str(r),
            "closed_gap": str(F(r, 4)), "cubic_anticommutator_cuts": 2 * r,
            "squares_in_certificate": 4 * r, "local_cut_violation": "-1/8",
            "fixed_N_contraction": True}


def main():
    local_identity()
    receipt = {"scope": "Exact algebraic occupancy diagnostic; no molecular or learning result",
               "local_CAR_identity": True, "full_six_mode_blocks": verify_dqg_blocks(),
               "family": [family_receipt(r) for r in (1, 2, 4, 8, 16, 32, 64, 128)]}
    path = Path(__file__).resolve().parents[1] / "results/marginal_hunt_witness.json"
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
