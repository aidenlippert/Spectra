"""Re-optimize the number ideal after fixing a truncated SOS factor budget.

The factors are treated as immutable.  A linear-programming pass chooses the
number-conserving multiplier and scalar ``b`` against the exact CAR
coefficient map, minimizing ``eta - b`` where ``eta`` is the coefficient-l1
residual.  The emitted certificate is then checked by the independent exact
verifier; floating point output is never used as a proof.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
from scipy import sparse
from scipy.optimize import linprog

from experiments.marginal_symbolic import (
    add, decode, encode, multiplier_basis, mono, number_shift, product,
    scale, verify,
)
from experiments.marginal_symbolic import expand_squares


def _map(poly, words):
    return np.array([float(poly.get(w, 0)) for w in words], dtype=float)


def reoptimize(source: dict, rounding: int = 10**9) -> tuple[dict, dict]:
    # Require a verifier-compatible source before optimizing.  Factor layout
    # is part of the certificate semantics and must never be guessed here.
    verify(source)
    source = dict(source)
    raw_blocks = source.get("blocks", [])
    source["blocks"] = []
    for block in raw_blocks:
        b = dict(block)
        if b["factor"] and len(b["factor"][0]) != len(b["words"]):
            raise ValueError(f"factor shape {len(b['factor'])}x{len(b['factor'][0])} does not match word count {len(b['words'])} in block {b.get('name')}")
        source["blocks"].append(b)
    modes, particles = source["modes"], source["particles"]
    degree = source.get("operator_degree", 3)
    h = decode(source["hamiltonian"], modes, 4)
    squares, stats = expand_squares(source["blocks"], source["denominator"], modes, degree)
    target = add(h, scale(squares, -1))
    basis = multiplier_basis(modes, max_body=degree - 1)
    shift = number_shift(modes, particles)
    ideal_cols = [product(shift, q) for q in basis]
    words = sorted(set(target) | {w for q in ideal_cols for w in q}, key=lambda w: (len(w), w))
    rhs = _map(target, words)
    A = np.column_stack([_map(q, words) for q in ideal_cols] + [_map(mono(()), words)])
    # residual = rhs - A @ [x,b].  t >= +/- residual; objective sum(t)-b.
    n = A.shape[1]
    Aub = sparse.vstack([
        sparse.hstack([sparse.csr_matrix(A), -sparse.eye(len(words))]),
        sparse.hstack([-sparse.csr_matrix(A), -sparse.eye(len(words))]),
    ], format="csr")
    bub = np.concatenate([rhs, -rhs])
    c = np.concatenate([np.zeros(n - 1), [-1.0], np.ones(len(words))])
    result = linprog(c, A_ub=Aub, b_ub=bub,
                     bounds=[(None, None)] * n + [(0, None)] * len(words),
                     method="highs")
    if not result.success:
        raise RuntimeError(f"ideal LP failed: {result.message}")
    den = int(rounding)
    # The first n entries are [multiplier coefficients..., b].  The remaining
    # entries are LP absolute-value slacks and are not certificate data.
    coeff = [F(int(round(v * den)), den) for v in result.x[:n]]
    x = add(*(scale(q, a) for q, a in zip(basis, coeff[:-1]) if a))
    certificate = dict(source)
    certificate["number_multiplier"] = encode(x)
    certificate.pop("multiplier_orbits", None)
    certificate["b"] = str(coeff[-1])
    receipt = verify(certificate)
    receipt.update({
        "method": "fixed_factors_l1_ideal_repair",
        "lp_status": int(result.status),
        "lp_objective_float": float(result.fun),
        "lp_iterations": int(result.nit if result.nit is not None else 0),
        "multiplier_basis_dimension": len(basis),
        "coefficient_rows": len(words),
        "multiplier_rounding_denominator": den,
        "fixed_factor_nonzeros": stats["factor_nonzeros"],
        "certificate_bytes": len(json.dumps(certificate, separators=(",", ":")).encode()),
    })
    return certificate, receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--rounding", type=int, default=10**9)
    args = parser.parse_args()
    started = time.monotonic()
    source = json.loads(args.input.read_text())
    certificate, receipt = reoptimize(source, args.rounding)
    receipt["wall_seconds"] = time.monotonic() - started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
