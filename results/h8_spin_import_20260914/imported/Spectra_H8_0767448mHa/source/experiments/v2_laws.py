"""Finite-shot discovery of reusable binary linear response laws over GF(2)."""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import comb
from hashlib import sha256
from math import ceil, log, isfinite
from typing import Sequence

@dataclass(frozen=True)
class CalibrationRecord:
    context: tuple[int, ...]
    outcomes: tuple[int, ...]

@dataclass(frozen=True)
class LearningResult:
    mask: tuple[int, ...] | None
    rank: int
    certificate_valid: bool
    confidence_delta: float
    shots_per_context: int
    xor_operations: int
    data_hash: str
    exact_error_bound: str = "0"
    reason: str = ""

def required_shots(d: int, eta: float, delta: float) -> int:
    _validate_params(d, eta, delta)
    if eta == 0: return 1
    k = ceil(log(d / delta) / (2 * (0.5 - eta) ** 2))
    k = k if k % 2 else k + 1
    while _exact_union_bound(d, k, eta) > _frac(delta):
        k += 2
    return k

def _validate_params(d, eta, delta):
    if type(d) is not int or not 1 <= d <= 256: raise ValueError("d must be positive")
    if isinstance(eta, bool) or not isfinite(eta) or not 0 <= eta < 0.5: raise ValueError("eta must satisfy 0 <= eta < 1/2")
    if isinstance(delta, bool) or not isfinite(delta) or not 0 < delta < 1: raise ValueError("delta must satisfy 0 < delta < 1")

def _bit(v):
    if isinstance(v, bool) or not isinstance(v, int) or v not in (0, 1): raise ValueError("bits must be integers 0 or 1")

def _frac(v):
    return v if isinstance(v, Fraction) else Fraction(str(v))

def _exact_union_bound(records_count, k, eta):
    e = _frac(eta); q = Fraction(0)
    for j in range((k + 1)//2, k + 1): q += comb(k, j) * e**j * (1-e)**(k-j)
    return records_count * q

def _rank_solve(rows: list[tuple[list[int], int]]) -> tuple[int, tuple[int, ...] | None, int, bool]:
    d = len(rows[0][0]); r = 0; ops = 0
    for col in range(d):
        pivot = next((i for i in range(r, len(rows)) if rows[i][0][col]), None)
        if pivot is None: continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][0][col]:
                rows[i] = ([a ^ b for a, b in zip(rows[i][0], rows[r][0])], rows[i][1] ^ rows[r][1]); ops += d + 1
        r += 1
    for row, rhs in rows:
        if not any(row) and rhs: return r, None, ops, False
    if r < d: return r, None, ops, True
    solution = [0] * d
    for row, rhs in rows:
        if not any(row): continue
        pivot = next(i for i, x in enumerate(row) if x)
        solution[pivot] = rhs
    return r, tuple(solution), ops, True

def learn(records: Sequence[CalibrationRecord], eta: float, delta: float) -> LearningResult:
    if not records: raise ValueError("records must be nonempty")
    d = len(records[0].context); k = len(records[0].outcomes)
    _validate_params(d, eta, delta)
    if k < 1 or k % 2 == 0: raise ValueError("outcomes require a positive odd shot count")
    rows = []
    for rec in records:
        if len(rec.context) != d or len(rec.outcomes) != k: raise ValueError("uniform record lengths required")
        for v in rec.context + rec.outcomes: _bit(v)
        majority = 1 if sum(rec.outcomes) > k // 2 else 0
        rows.append((list(rec.context), majority))
    rank, mask, ops, consistent = _rank_solve(rows)
    payload = repr(tuple((r.context, r.outcomes) for r in records)).encode()
    digest = sha256(payload).hexdigest()
    exact = _exact_union_bound(len(records), k, eta)
    bound = float(exact)
    valid = mask is not None and consistent and exact <= _frac(delta)
    reason = "inconsistent majority constraints" if not consistent else ("rank-deficient calibration" if mask is None else "")
    return LearningResult(mask, rank, valid, bound, k, ops, digest, str(exact), reason)

def verify_certificate(records: Sequence[CalibrationRecord], result: LearningResult, eta: float, delta: float) -> bool:
    """Check model-conditional procedure coverage, not posterior certainty.

    The algebra checker does not call learn. The same exact binomial-tail helper
    is used, and is independently checked against exhaustive short bit strings
    in tests. Physical source assumptions cannot be certified from these records.
    """
    try:
        if not records or result.mask is None: return False
        d = len(records[0].context); k = len(records[0].outcomes)
        _validate_params(d, eta, delta)
        if k < 1 or k % 2 == 0 or type(result.shots_per_context) is not int or result.shots_per_context != k:
            return False
        if type(result.mask) is not tuple or len(result.mask) != d:
            return False
        for v in result.mask: _bit(v)
        if result.data_hash != sha256(repr(tuple((r.context, r.outcomes) for r in records)).encode()).hexdigest():
            return False
        rows = []
        for rec in records:
            if len(rec.context) != d or len(rec.outcomes) != k: return False
            for v in rec.context + rec.outcomes: _bit(v)
            label = int(sum(rec.outcomes) > k // 2)
            rows.append((list(rec.context), label))
            if predict(result.mask, rec.context) != label: return False
        rank, _, ops, consistent = _rank_solve(rows)
        exact = _exact_union_bound(len(records), k, eta)
        return (consistent and rank == d and type(result.rank) is int and result.rank == d
                and result.certificate_valid is True and result.reason == ""
                and type(result.xor_operations) is int and result.xor_operations == ops
                and result.confidence_delta == float(exact)
                and result.exact_error_bound == str(exact) and exact <= _frac(delta))
    except (AttributeError, TypeError, ValueError, OverflowError):
        return False

def predict(mask: Sequence[int], inputs: Sequence[int]) -> int:
    if len(mask) != len(inputs): raise ValueError("mask and inputs must have equal length")
    for v in tuple(mask) + tuple(inputs): _bit(v)
    return sum(a & b for a, b in zip(mask, inputs)) % 2
