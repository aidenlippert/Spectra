"""Finite-patch CSS stabilizer support diagnostic over GF(2).

This checks a necessary local-generation condition related to TQO-2.  It is
deliberately a finite linear-algebra test; it does not prove TQO-2 or thermal
memory lifetime.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence


def _bits(row: Iterable[int], n: int) -> tuple[int, ...]:
    values = tuple(row)
    if len(values) != n or any(type(x) is not int or x not in (0, 1) for x in values):
        raise ValueError("rows must be binary vectors of the declared length")
    return values


def _rref(rows: Sequence[Sequence[int]], n: int) -> tuple[list[list[int]], list[int]]:
    a = [list(_bits(row, n)) for row in rows]
    pivots: list[int] = []
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        for i in range(len(a)):
            if i != r and a[i][c]:
                a[i] = [x ^ y for x, y in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == len(a):
            break
    return a[:r], pivots


def _nullspace(rows: Sequence[Sequence[int]], n: int) -> list[tuple[int, ...]]:
    rref, pivots = _rref(rows, n)
    free = [c for c in range(n) if c not in pivots]
    result: list[tuple[int, ...]] = []
    for f in free:
        v = [0] * n
        v[f] = 1
        for i, p in enumerate(pivots):
            v[p] = rref[i][f]
        result.append(tuple(v))
    return result


def _span_intersection_supported(rows: Sequence[Sequence[int]], allowed: set[int], n: int) -> list[tuple[int, ...]]:
    """Basis for span(rows) vectors whose support is contained in allowed."""
    outside = [c for c in range(n) if c not in allowed]
    constraints = [[row[c] for row in rows] for c in outside]
    coefficients = _nullspace(constraints, len(rows))
    vectors = []
    for coeff in coefficients:
        vectors.append(tuple(sum(coeff[i] * rows[i][c] for i in range(len(rows))) % 2 for c in range(n)))
    basis, _ = _rref(vectors, n)
    return [tuple(row) for row in basis]


def _in_span(vector: tuple[int, ...], rows: Sequence[Sequence[int]], n: int) -> bool:
    before, pivots = _rref(rows, n)
    rank_before = len(pivots)
    after, pivots_after = _rref([*rows, vector], n)
    return len(pivots_after) == rank_before


@dataclass(frozen=True)
class Diagnostic:
    passed: bool
    witness: tuple[int, ...] | None = None
    sector: str | None = None


def _check_sector(generators: Sequence[Sequence[int]], region_a: set[int], buffer_b: set[int], n: int, sector: str) -> Diagnostic:
    global_rows = [_bits(row, n) for row in generators]
    local_rows = [row for row in global_rows if all(row[c] == 0 for c in range(n) if c not in buffer_b)]
    global_a = _span_intersection_supported(global_rows, region_a, n)
    local_a = _span_intersection_supported(local_rows, region_a, n)
    for witness in global_a:
        if any(witness) and not _in_span(witness, local_a, n):
            return Diagnostic(False, witness, sector)
    return Diagnostic(True)


def diagnose_css(
    x_generators: Sequence[Sequence[int]],
    z_generators: Sequence[Sequence[int]],
    n_qubits: int,
    region_a: Iterable[int],
    buffer_b: Iterable[int],
) -> tuple[Diagnostic, Diagnostic]:
    """Compare global-vs-buffer-supported X and Z stabilizer support spaces."""
    if type(n_qubits) is not int or n_qubits <= 0:
        raise ValueError("n_qubits must be positive")
    a, b = set(region_a), set(buffer_b)
    if any(type(q) is not int for q in a | b) or not a <= b or any(q < 0 or q >= n_qubits for q in b):
        raise ValueError("require region_a subset buffer_b within the patch")
    x = [_bits(row, n_qubits) for row in x_generators]
    z = [_bits(row, n_qubits) for row in z_generators]
    if any(sum(xr[q] * zr[q] for q in range(n_qubits)) % 2 for xr in x for zr in z):
        raise ValueError("CSS X and Z generators must commute over GF(2)")
    return (_check_sector(x, a, b, n_qubits, "X"), _check_sector(z, a, b, n_qubits, "Z"))
