"""Exact coefficient witnesses distinguish sector nulls from truncated-ideal equivalence.

The free scalar energy coordinate is deliberately not an ideal generator.
"""
from fractions import Fraction as F

from experiments.marginal_hunt_car import add, adj, mono, scale
from experiments.marginal_symbolic import canonical, number_shift, product


class IdealSpan:
    def __init__(self, columns):
        self.columns = []
        self.pivots = {}
        for label, polynomial in columns:
            if label in ("energy", "free_scalar", "identity"):
                raise ValueError("The energy coordinate is not a vanishing ideal")
            p = canonical(polynomial)
            if p and set(p) == {()}:
                raise ValueError("A nonzero scalar cannot be a vanishing ideal")
            index = len(self.columns)
            self.columns.append((label, p))
            remainder, witness = self.reduce(p)
            if remainder:
                pivot = min(remainder, key=lambda w: (len(w), w))
                v = remainder[pivot]
                # remainder = column[index] - sum(witness[i] column[i])
                coordinate = {index: F(1, 1) / v, **{i: -c / v for i, c in witness.items()}}
                self.pivots[pivot] = (scale(remainder, 1 / v), coordinate)
        if not self.reduce(mono(()))[0]:
            raise ValueError("Declared vanishing columns span the energy scalar")

    def reduce(self, polynomial):
        remainder = dict(canonical(polynomial))
        witness = {}
        # A pivot row has no earlier pivot coefficients, so one ordered pass suffices.
        for pivot in sorted(self.pivots, key=lambda w: (len(w), w)):
            c = remainder.get(pivot, 0)
            if not c:
                continue
            row, coordinates = self.pivots[pivot]
            remainder = add(remainder, scale(row, -c))
            for i, a in coordinates.items():
                witness[i] = witness.get(i, 0) + c * a
        return remainder, {i: c for i, c in witness.items() if c}

    def witness(self, polynomial):
        p = canonical(polynomial)
        remainder, coordinates = self.reduce(p)
        reconstructed = add(*(scale(self.columns[i][1], c) for i, c in coordinates.items()))
        if add(p, scale(reconstructed, -1)) != remainder:
            raise AssertionError("Broken exact coefficient witness")
        return {"equivalent_in_declared_ideal": not remainder,
                "coordinates": coordinates, "unresolved_polynomial": remainder}


def number_null_cross_terms(modes, particles, reduced, multiplier, span):
    """B = reduced + multiplier*(N-n); every Gram cross term is checked.

    Physical equality on the fixed-N input sector is unconditional. Equality
    under a particular truncated optimization is a separate exact span test.
    """
    constraint = number_shift(modes, particles)
    null = product(multiplier, constraint)
    original = add(reduced, null)
    diagonal = add(product(canonical(adj(original)), original),
                   scale(product(canonical(adj(reduced)), reduced), -1))
    # A Hermitian square difference includes both cross terms and null^dagger null.
    result = span.witness(diagonal)
    return {"physical_fixed_N_equality": True, "includes_both_cross_terms_and_null_square": True,
            "original": original, "reduced": reduced, "difference": diagonal, **result,
            "optimization_status": "equivalent" if result["equivalent_in_declared_ideal"]
            else "physical_strengthening_not_current_optimization_equivalence"}


def verify_witness(polynomial, columns, coordinates):
    if any(type(i) is not int or not 0 <= i < len(columns) for i in coordinates):
        return False
    if any(label in ("energy", "identity", "free_scalar") for label, _ in columns):
        return False
    if any(p and set(canonical(p)) == {()} for _, p in columns):
        return False
    # The independent reconstruction entry point must also reject combinations
    # of apparently non-scalar columns that secretly make the energy free.
    try:
        IdealSpan(columns)
    except ValueError:
        return False
    return canonical(polynomial) == add(*(scale(columns[i][1], F(c)) for i, c in coordinates.items()))
