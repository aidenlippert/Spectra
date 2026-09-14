"""Small exact Pauli-string closure utilities for intervention tests."""
from __future__ import annotations

from fractions import Fraction
from typing import Mapping, Iterable

_MUL = {
    ("I", "I"): (1, "I"), ("I", "X"): (1, "X"), ("I", "Y"): (1, "Y"), ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"), ("Y", "I"): (1, "Y"), ("Z", "I"): (1, "Z"),
    ("X", "X"): (1, "I"), ("Y", "Y"): (1, "I"), ("Z", "Z"): (1, "I"),
    ("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"),
    ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
    ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y"),
}

def _check_label(label: str) -> None:
    if not isinstance(label, str) or not label or any(c not in "IXYZ" for c in label):
        raise ValueError("Pauli labels must be nonempty strings over I,X,Y,Z")

def multiply(a: str, b: str) -> tuple[complex, str]:
    _check_label(a); _check_label(b)
    if len(a) != len(b): raise ValueError("Pauli strings must have equal length")
    phase, out = 1, []
    for x, y in zip(a, b):
        p, z = _MUL[(x, y)]; phase *= p; out.append(z)
    return phase, "".join(out)

def commutator_i(hamiltonian: Mapping[str, Fraction | int], observable: str) -> dict[str, Fraction]:
    """Return exact real Pauli coefficients of i[H,O] for Hermitian H."""
    _check_label(observable); n = len(observable); out: dict[str, Fraction] = {}
    for label, coeff in hamiltonian.items():
        _check_label(label)
        if len(label) != n: raise ValueError("all Pauli strings must have equal length")
        c = Fraction(coeff)
        p1, z1 = multiply(label, observable); p2, z2 = multiply(observable, label)
        if z1 != z2: raise AssertionError("Pauli products disagree")
        phase_difference = 1j * (p1 - p2)
        if abs(phase_difference.imag) > 0:
            raise ValueError("Hermitian commutator produced a non-real coefficient")
        value = c * int(phase_difference.real)
        if value:
            out[z1] = out.get(z1, Fraction(0)) + Fraction(value)
    return {k: v for k, v in out.items() if v}

def closure(controls: Iterable[Mapping[str, Fraction | int]], seeds: Iterable[str], budget: int = 128) -> tuple[set[str], bool]:
    controls = list(controls); basis = set(seeds)
    if budget < 1: return set(), False
    for s in basis: _check_label(s)
    n = len(next(iter(basis))) if basis else None
    if any(len(s) != n for s in basis): raise ValueError("uniform Pauli lengths required")
    for h in controls:
        _check_basis(list(h))
        if n is not None and any(len(s) != n for s in h): raise ValueError("uniform Pauli lengths required")
    if len(basis) > budget: return set(sorted(basis)[:budget]), False
    changed = True
    while changed:
        changed = False
        for h in controls:
            for o in sorted(basis):
                for z in sorted(commutator_i(h, o)):
                    if z not in basis:
                        if len(basis) >= budget: return basis, False
                        basis.add(z); changed = True
    return basis, True

def generator_matrix(hamiltonian: Mapping[str, Fraction | int], basis: Iterable[str]) -> tuple[list[str], list[list[Fraction]], dict[str, Fraction]]:
    """Return row-oriented projected generator: zdot[i] = sum_j A[i][j] z[j]."""
    labels = list(dict.fromkeys(basis)); _check_basis(labels)
    index = {x:i for i,x in enumerate(labels)}; matrix = [[Fraction(0) for _ in labels] for _ in labels]; residual = {}
    for j, o in enumerate(labels):
        deriv = commutator_i(hamiltonian, o)
        residual[o] = sum(abs(v) for z,v in deriv.items() if z not in index)
        for z,v in deriv.items():
            if z in index: matrix[j][index[z]] = v
    return labels, matrix, residual

def _check_basis(labels: list[str]) -> None:
    if not labels: raise ValueError("basis must be nonempty")
    for x in labels: _check_label(x)
    if len({len(x) for x in labels}) != 1: raise ValueError("uniform Pauli lengths required")

def commutation_adjacency(terms: Iterable[str]) -> int:
    labels = list(terms); _check_basis(labels); return sum(1 for i,a in enumerate(labels) for b in labels[i+1:] if multiply(a,b)[0] != multiply(b,a)[0])
