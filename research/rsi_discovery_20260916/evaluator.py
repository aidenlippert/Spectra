"""Frozen acceptance authority. Candidate code never runs in this process.

The molecular gate reuses the existing exact Cholesky-remainder checker and
Rayleigh quotient on every block of the entire fixed-N sector. Enumeration is
explicit and capped. These are finite rational model certificates, not physical
model validation or a replacement for the large-system SOS/MPS pipeline.
"""
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import sys

from research.transfer_followup_20260915.enumerated_baseline import blocks, shifted, quotient
from research.compact_response_20260913.closure import check_factor
from research.rsi_discovery_20260916.ledger import digest

TARGET = F(1, 625)


def validate_fixture(data):
    m, n = data["modes"], data["particles"]
    if type(m) is not int or type(n) is not int or not 1 <= m <= 16 or not 0 <= n <= m:
        raise ValueError("Invalid or unsupported finite fermionic sector")
    if comb(m, n) > 2000:
        raise ValueError("Enumeration envelope exceeded")
    if not data.get("hamiltonian"):
        raise ValueError("Nonempty exact Hamiltonian required")
    return blocks(data)


def constructor_payload(data):
    groups, den, cost = validate_fixture(data)
    return {"blocks": [{"alpha": a, "H": h} for a, states, h in groups], "den": den}, cost


def verify_certificate(data, cert):
    if cert.get("fixture_sha256") != digest(data):
        raise ValueError("Exact task identity mismatch")
    if cert.get("kind") != "discovered_full_N_certificate_v1":
        raise ValueError("Unsupported certificate type")
    groups, den, cost = validate_fixture(data)
    rows = cert["blocks"]
    if len(rows) != len(groups) or [r["alpha_count"] for r in rows] != [a for a, s, h in groups]:
        raise ValueError("Missing, duplicated or reordered spin-projection block")
    lows, ups, margins = [], [], []
    for (a, states, h), row in zip(groups, rows):
        if type(row["alpha_count"]) is not int or not isinstance(row["lower_Ha"], str):
            raise ValueError("Malformed exact endpoint")
        lower = F(row["lower_Ha"])
        fd = row["factor_denominator"]
        if type(fd) is not int or not 0 < fd <= 10**16:
            raise ValueError("Positive bounded integer factor denominator required")
        if any(abs(x) > 10**24 for r in row["factor"] for x in r if type(x) is int):
            raise ValueError("Factor bit-size envelope exceeded")
        k, kd = shifted(h, den, lower)
        margin = check_factor(k, kd, F(0), row["factor"], fd)
        upper = quotient(h, den, row["upper_vector"])
        if lower > upper:
            raise ValueError("Lower endpoint exceeds concrete Rayleigh quotient")
        lows.append(lower)
        ups.append(upper)
        margins.append(str(margin))
    lower, upper = min(lows), min(ups)
    if upper - lower > TARGET:
        raise ValueError("Frozen 1.6 mHa target missed")
    return {"status": "verified_instance", "lower_Ha": str(lower), "upper_Ha": str(upper),
            "width_mHa": float(1000 * (upper - lower)), "margins_Ha": margins,
            "assumptions": ["Supplied rational Hamiltonian", "entire fixed particle-number sector",
                            "real Hermitian, number and alpha-count conserving coefficients"],
            "cost": cost, "integral_error_certified": False, "physical_model_error_certified": False}


def reference_gram(factor):
    if not factor or len(factor) > 500 or any(len(r) != i + 1 for i, r in enumerate(factor)):
        raise ValueError("Bounded complete triangular integer factor required")
    if any(type(v) is not int or v.bit_length() > 256 for r in factor for v in r):
        raise ValueError("Exact bounded integers required")
    n = len(factor)
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            value = sum(x * y for x, y in zip(factor[i], factor[j]))
            out[i][j] = out[j][i] = value
    return out


def verify_grams(payloads, results):
    if not isinstance(results, list) or len(results) != len(payloads):
        raise ValueError("Incomplete component outputs")
    for payload, result in zip(payloads, results):
        n = len(payload["factor"])
        if (not isinstance(result, list) or len(result) != n or
                any(not isinstance(r, list) or len(r) != n or any(type(x) is not int for x in r) for r in result)):
            raise ValueError("Exact complete integer Gram result required")
        if result != reference_gram(payload["factor"]):
            raise ValueError("Exact equivalence failed")
    return {"status": "verified_instance", "instances": len(payloads),
            "scope": "Exact equality on the recorded finite input suite; not a universal program proof"}


def main():
    mode, input_path, output_path = sys.argv[1:]
    data = json.loads(Path(input_path).read_text())
    if mode == "certificate":
        result = verify_certificate(data["fixture"], data["certificate"])
    elif mode == "grams":
        result = verify_grams(data["payloads"], data["outputs"])
    else:
        raise ValueError("Unknown independent evaluation")
    if any(name in sys.modules for name in ("numpy", "scipy", "pyscf", "flint")):
        raise RuntimeError("Independent gate imported a numerical/candidate backend")
    result["numerical_imports"] = []
    Path(output_path).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
