"""Bounded six-mode certificate discovery. Exact export, numerical proposals.

Uses installed NumPy/SciPy/CVXPY for proposals. Verification below uses integer
Fock matrices and rational arithmetic, not the SDP solver's status or objective.
This sector representation is exponential and is only a diagnostic.
"""
from fractions import Fraction
from itertools import combinations
import json
from numbers import Rational
from pathlib import Path
import time

import numpy as np
import scipy.sparse as sp
import cvxpy as cp

MODES, PARTICLES = 6, 3
SECTOR = tuple(s for s in range(1 << MODES) if s.bit_count() == PARTICLES)


def dagger(word):
    return tuple((1 - c, i) for c, i in reversed(word))


def word_matrix(word):
    charge = sum(2 * c - 1 for c, _ in word)
    outputs = [s for s in range(1 << MODES) if s.bit_count() == PARTICLES + charge]
    indices = {s: i for i, s in enumerate(outputs)}
    matrix = np.zeros((len(outputs), len(SECTOR)), dtype=np.int64)
    for column, state in enumerate(SECTOR):
        value = 1
        for creation, mode in reversed(word):
            if ((state >> mode) & 1) == creation:
                value = 0
                break
            value *= (-1) ** ((state & ((1 << mode) - 1)).bit_count())
            state ^= 1 << mode
        if value:
            matrix[indices[state], column] = value
    return matrix


def hamiltonian(t, asymmetric=False):
    t = Fraction(t)
    matrix = np.zeros((20, 20), dtype=object)
    for col, state in enumerate(SECTOR):
        for triple in ((0, 1, 2), (3, 4, 5)):
            n = sum((state >> i) & 1 for i in triple)
            matrix[col, col] += n * (n - 1) // 2
    edges = [(0, 3, Fraction(1)), (1, 4, Fraction(1)), (2, 5, Fraction(1))]
    if asymmetric:
        edges = [(0, 3, Fraction(1)), (1, 4, Fraction(7, 10)),
                 (2, 5, Fraction(13, 10)), (0, 4, Fraction(2, 5))]
    for i, j, weight in edges:
        hop = word_matrix(((1, i), (0, j)))
        matrix -= t * weight * (hop + hop.T)
    return matrix


def make_block(name, words, mix=None):
    words = tuple(tuple(tuple(letter) for letter in word) for word in words)
    raw = np.array([word_matrix(w) for w in words])
    if mix is None:
        mix = np.eye(len(words))
    return {"name": name, "words": words, "raw": raw,
            "mix": np.asarray(mix, dtype=float)}


def baseline_blocks():
    linear = [((0, i),) for i in range(6)]
    pairs = [((0, j), (0, i)) for i, j in combinations(range(6), 2)]
    number = [((1, i), (0, j)) for i in range(6) for j in range(6)]
    return [make_block("linear-", linear), make_block("linear+", [dagger(w) for w in linear]),
            make_block("pair-", pairs), make_block("pair+", [dagger(w) for w in pairs]),
            make_block("particle-hole", number)]


def seed_blocks():
    # Independent squares are sound in the fixed-sector identity. They are not
    # individually asserted to define two-body inequalities.
    result = []
    for k, triple in enumerate(((0, 1, 2), (3, 4, 5))):
        word = tuple((0, i) for i in triple)
        result += [make_block(f"triple{k}-", [word]), make_block(f"triple{k}+", [dagger(word)])]
    return result


def cubic_pools(kind="density"):
    # Restrict to density-assisted annihilators; exclude the complete pure
    # triple dictionary that trivially spans the N=3 sector.
    words = [((0, i),) for i in range(6)]
    if kind == "density":
        words += [((1, j), (0, j), (0, i)) for i in range(6) for j in range(6) if i != j]
    elif kind == "mixed":
        words += [((1, k), (0, j), (0, i)) for i, j in combinations(range(6), 2) for k in range(6)]
    else:
        raise ValueError("Unknown dictionary")
    return [make_block(f"{kind}-", words), make_block(f"{kind}+", [dagger(w) for w in words])]


def solve(h, blocks):
    b = cp.Variable()
    expression = b * np.eye(20)
    variables = []
    shared = {}
    for block in blocks:
        w = np.einsum("uv,vai->uai", block["mix"], block["raw"])
        k = len(w)
        tie = block.get("tie")
        if tie is not None and tie in shared:
            variable = shared[tie]
            if variable.shape != (k, k):
                raise ValueError("Tied Gram blocks must have equal dimensions")
        else:
            variable = cp.Variable((k, k), PSD=True)
            if tie is not None:
                shared[tie] = variable
        products = np.einsum("uai,vaj->ijuv", w, w).reshape(400, k * k)
        expression = expression + cp.reshape(sp.csr_matrix(products) @ cp.reshape(variable, (k*k,), order="C"), (20, 20), order="C")
        variables.append(variable)
    match = expression == np.asarray(h, dtype=float)
    problem = cp.Problem(cp.Maximize(b), [match])
    start = time.monotonic()
    problem.solve(solver="SCS", eps=2e-7, max_iters=20000)
    elapsed = time.monotonic() - start
    if b.value is None or any(q.value is None for q in variables):
        raise RuntimeError(f"No certificate proposal: {problem.status}")
    y = np.asarray(match.dual_value)
    if np.trace(y) < 0:
        y = -y
    return {"b": float(b.value), "grams": [q.value for q in variables],
            "pseudo_state": (y + y.T) / 2, "seconds": elapsed, "status": problem.status}


def export_certificate(h, blocks, solution, digits=6):
    denominator = 10 ** digits
    factors = []
    for block, gram in zip(blocks, solution["grams"]):
        values, vectors = np.linalg.eigh((gram + gram.T) / 2)
        # Clipping affects tightness only: the exact factor, not floating Q,
        # defines the accepted PSD operator.
        l = (np.sqrt(np.maximum(values, 0))[:, None] * vectors.T) @ block["mix"]
        ints = [[int(round(x * denominator)) for x in row] for row in l]
        ints = [row for row in ints if any(row)]
        factors.append({"name": block["name"], "words": block["words"], "factor": ints})
    vals, vecs = np.linalg.eigh(np.asarray(h, dtype=float))
    trial = [int(round(x * denominator)) for x in vecs[:, 0]]
    certificate = {"mode_count": 6, "particle_number": 3, "denominator": denominator,
                   "b_numerator": int(round(solution["b"] * denominator ** 2)),
                   "blocks": factors, "trial": trial}
    return certificate, float(vals[0])


def verify_certificate(h, certificate):
    """Exact finite-sector replay. Returns bounds even for a poor proposal."""
    if h.shape != (20, 20) or any(not isinstance(x, Rational) for x in h.flat):
        raise ValueError("An exact rational sector Hamiltonian is required")
    if not np.array_equal(h, h.T):
        raise ValueError("Non-Hermitian Hamiltonian")
    if certificate["mode_count"] != 6 or certificate["particle_number"] != 3:
        raise ValueError("Unsupported sector")
    den = certificate["denominator"]
    if type(den) is not int or den <= 0:
        raise ValueError("Invalid denominator")
    total = np.zeros((20, 20), dtype=object)
    for block in certificate["blocks"]:
        words = block["words"]
        raw = []
        charges = set()
        for word in words:
            if not word or len(word) > 3 or any(type(c) is not int or c not in (0, 1) or type(m) is not int or m not in range(6) for c, m in word):
                raise ValueError("Invalid operator word")
            charges.add(sum(2*c-1 for c, _ in word))
            raw.append(word_matrix(word).astype(object))
        if len(charges) != 1:
            raise ValueError("Mixed charges in Gram block")
        for row in block["factor"]:
            if len(row) != len(words) or any(type(v) is not int for v in row):
                raise ValueError("Invalid rational factor")
            f = sum((v * w for v, w in zip(row, raw)), np.zeros_like(raw[0]))
            total += f.T @ f
    bnum = certificate["b_numerator"]
    if type(bnum) is not int:
        raise ValueError("Invalid lower-bound scalar")
    residual = h - (total + bnum * np.eye(20, dtype=object)) / Fraction(den**2)
    if not np.array_equal(residual, residual.T):
        raise ValueError("Non-Hermitian residual")
    eta = max(sum(abs(x) for x in row) for row in residual)
    lower = Fraction(bnum, den**2) - eta
    trial = certificate["trial"]
    if len(trial) != 20 or any(type(v) is not int for v in trial) or not any(trial):
        raise ValueError("Invalid variational state")
    x = np.asarray(trial, dtype=object)
    upper = Fraction(x @ h @ x, x @ x)
    return {"lower": str(lower), "upper": str(upper), "width": str(upper-lower),
            "residual_row_norm": str(eta), "lower_float": float(lower),
            "upper_float": float(upper), "width_float": float(upper-lower),
            "factor_rows": sum(len(z["factor"]) for z in certificate["blocks"])}


def adaptive(h, max_rounds=6, pool_kind="density", paired=False):
    blocks = baseline_blocks() + seed_blocks()
    pools = cubic_pools(pool_kind)
    selected = [[] for _ in pools]
    history = []
    for iteration in range(max_rounds + 1):
        current = blocks + [make_block(p["name"], p["words"], rows) for p, rows in zip(pools, selected) if rows]
        solution = solve(h, current)
        history.append({"iteration": iteration, "numeric_lower": solution["b"],
                        "selected_directions": sum(map(len, selected)), "seconds": solution["seconds"]})
        if iteration == max_rounds:
            break
        improved = False
        moments = []
        for k, p in enumerate(pools):
            w = p["raw"]
            moment = np.einsum("uai,ij,vaj->uv", w, solution["pseudo_state"], w)
            moments.append((moment+moment.T)/2)
        if paired:
            moment = moments[0] + moments[1]
            eig, vectors = np.linalg.eigh(moment)
            history[-1]["minimum_paired_moment"] = float(eig[0])
            if eig[0] < -1e-7:
                for rows in selected:
                    rows.append(vectors[:, 0])
                improved = True
        else:
          for k, moment in enumerate(moments):
            eig, vectors = np.linalg.eigh((moment+moment.T)/2)
            if eig[0] < -1e-7:
                selected[k].append(vectors[:, 0])
                improved = True
        if not improved:
            break
    return current, solution, history


def main():
    cases = [("matched", Fraction(1, 5), False), ("matched", Fraction(1), False),
             ("asymmetric", Fraction(1, 5), True), ("asymmetric", Fraction(1), True)]
    receipts = []
    out = Path(__file__).resolve().parents[1] / "results/marginal_hopping"
    out.mkdir(exist_ok=True)
    for label, t, asym in cases:
        h = hamiltonian(t, asym)
        for arm in ("quadratic", "local_cubic", "adaptive", "density_pool"):
            if arm == "adaptive":
                blocks, proposal, history = adaptive(h)
            else:
                blocks = baseline_blocks()
                if arm != "quadratic":
                    blocks += seed_blocks()
                if arm == "density_pool":
                    blocks += cubic_pools()
                proposal = solve(h, blocks)
                history = []
            certificate, reference = export_certificate(h, blocks, proposal)
            certificate.update({"t": str(t), "asymmetric": asym})
            receipt = verify_certificate(h, certificate)
            receipt.update({"case": label, "t": str(t), "arm": arm, "reference_numeric": reference,
                            "proposal_seconds": proposal["seconds"], "history": history})
            filename = f"{label}_{str(t).replace('/', '_')}_{arm}.json"
            (out / filename).write_text(json.dumps(certificate) + "\n")
            receipt["certificate_file"] = filename
            receipts.append(receipt)
            print(json.dumps(receipt), flush=True)
            (out / "summary.json").write_text(json.dumps(receipts, indent=2) + "\n")


def recycled_subspace(h, rounds=3):
    base, pools = baseline_blocks()+seed_blocks(), cubic_pools("mixed")
    mix, history = np.empty((0, 96)), []
    for iteration in range(rounds+1):
        blocks = base + [make_block(p["name"], p["words"], mix) for p in pools] if len(mix) else base
        proposal = solve(h, blocks)
        cert, _ = export_certificate(h, blocks, proposal)
        result = verify_certificate(h, cert)
        history.append({"iteration": iteration, "directions_per_charge": len(mix),
                        "width": result["width_float"], "numeric_lower": proposal["b"],
                        "seconds": proposal["seconds"]})
        if iteration == rounds:
            break
        moment = sum(np.einsum("uai,ij,vaj->uv", p["raw"], proposal["pseudo_state"], p["raw"]) for p in pools)
        eig, vectors = np.linalg.eigh((moment+moment.T)/2)
        extra = vectors[:, eig < -1e-7].T
        if not len(extra):
            break
        _, values, vh = np.linalg.svd(np.vstack([mix, extra]), full_matrices=False)
        mix = vh[values > 1e-7]
    return blocks, proposal, history


def run_stage(stage):
    """Reproducible entry points for the additional bounded search controls."""
    out = Path(__file__).resolve().parents[1] / "results/marginal_hopping"
    out.mkdir(exist_ok=True)
    if stage == "replay":
        replay = []
        for path in sorted(out.glob("*.json")):
            data = json.loads(path.read_text())
            if not isinstance(data, dict) or "blocks" not in data:
                continue
            result = verify_certificate(hamiltonian(Fraction(data["t"]), data["asymmetric"]), data)
            if Fraction(result["width"]) < 0:
                raise AssertionError("Inconsistent certified interval")
            replay.append({"file": path.name, **result})
        (out / "replay.json").write_text(json.dumps(replay, indent=2) + "\n")
        print(f"Replayed {len(replay)} rational certificates")
        return
    records = []
    cases = [("matched", Fraction(1, 5), False), ("matched", Fraction(1), False),
             ("asymmetric", Fraction(1, 5), True), ("asymmetric", Fraction(1), True)]
    for label, t, asym in cases:
        if stage in ("subspace", "full", "anticommutator", "recycled") and t != Fraction(1, 5):
            continue
        h = hamiltonian(t, asym)
        base, pools = baseline_blocks() + seed_blocks(), cubic_pools("mixed")
        experiments = []
        if stage == "paired":
            blocks, proposal, history = adaptive(h, max_rounds=8, pool_kind="mixed", paired=True)
            experiments.append(("paired_mixed", blocks, proposal, {"history": history}))
        elif stage == "random":
            mix = np.random.default_rng(1729).normal(size=(8, 96))
            mix /= np.linalg.norm(mix, axis=1)[:, None]
            blocks = base + [make_block(p["name"], p["words"], mix) for p in pools]
            experiments.append(("random", blocks, solve(h, blocks), {"seed": 1729}))
        elif stage in ("full", "anticommutator"):
            if stage == "anticommutator":
                for p in pools:
                    p["tie"] = "anticommutator"
            blocks = base + pools
            experiments.append(("full_mixed" if stage == "full" else "full_anticommutator", blocks, solve(h, blocks), {}))
        elif stage == "recycled":
            blocks, proposal, history = recycled_subspace(h)
            experiments.append(("recycled_subspace", blocks, proposal, {"history": history}))
        elif stage == "subspace":
            initial = solve(h, base)
            moment = sum(np.einsum("uai,ij,vaj->uv", p["raw"], initial["pseudo_state"], p["raw"]) for p in pools)
            eig, vectors = np.linalg.eigh((moment + moment.T)/2)
            negative = eig < -1e-7
            chosen = vectors[:, negative].T
            for arm in ("negative_subspace", "random_equal_subspace"):
                mix = chosen.copy()
                if arm.startswith("random"):
                    mix = np.random.default_rng(1729).normal(size=chosen.shape)
                    mix /= np.linalg.norm(mix, axis=1)[:, None]
                blocks = base + [make_block(p["name"], p["words"], mix) for p in pools] if len(mix) else base
                experiments.append((arm, blocks, solve(h, blocks),
                                    {"directions_per_charge": len(mix), "negative_eigenvalues": eig[negative].tolist()}))
        else:
            raise ValueError("Unknown stage")
        for arm, blocks, proposal, details in experiments:
            cert, reference = export_certificate(h, blocks, proposal)
            cert.update({"t": str(t), "asymmetric": asym})
            result = verify_certificate(h, cert)
            filename = f"{label}_{str(t).replace('/', '_')}_{arm}.json"
            (out / filename).write_text(json.dumps(cert) + "\n")
            result.update({"case": label, "t": str(t), "arm": arm,
                           "reference_numeric": reference, "status": proposal["status"],
                           "proposal_seconds": proposal["seconds"], "certificate_file": filename, **details})
            records.append(result)
            (out / f"{stage}_summary.json").write_text(json.dumps(records, indent=2) + "\n")
            print(json.dumps(result), flush=True)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("baseline", "paired", "random", "full", "subspace", "recycled", "anticommutator", "replay"), default="baseline")
    stage = parser.parse_args().stage
    main() if stage == "baseline" else run_stage(stage)
