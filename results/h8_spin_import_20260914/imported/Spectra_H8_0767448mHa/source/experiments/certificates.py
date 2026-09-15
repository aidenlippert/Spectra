"""Exact rational anticommuting-group energy bounds and a separate checker."""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from typing import Iterable, Sequence


@dataclass(frozen=True)
class PauliTerm:
    pauli: str
    coeff: Fraction

    def __post_init__(self):
        if not isinstance(self.pauli, str) or not self.pauli or any(c not in 'IXYZ' for c in self.pauli):
            raise ValueError('invalid Pauli string')
        object.__setattr__(self, 'coeff', Fraction(self.coeff))


def _anti(a, b):
    if len(a) != len(b):
        raise ValueError('different Pauli widths')
    return sum(x != 'I' and y != 'I' and x != y for x, y in zip(a, b)) % 2 == 1


def _sqrt_interval(value, digits=12):
    if not isinstance(digits, int) or not 0 <= digits <= 100 or value < 0:
        raise ValueError('invalid precision or radicand')
    scale = 10 ** digits
    numerator = value.numerator * scale ** 2
    denominator = value.denominator
    lower = isqrt(numerator // denominator)
    upper = lower if lower ** 2 * denominator == numerator else lower + 1
    return Fraction(lower, scale), Fraction(upper, scale)


def _validated(terms, state):
    terms = list(terms)
    if not terms:
        raise ValueError('empty Hamiltonian')
    width = len(terms[0].pauli)
    if any(len(t.pauli) != width for t in terms):
        raise ValueError('different Pauli widths')
    state = tuple((0,) * width if state is None else state)
    if len(state) != width or any(type(x) is not int or x not in (0, 1) for x in state):
        raise ValueError('invalid computational basis state')
    return terms, state


def _expectation(term, state):
    if any(c in 'XY' for c in term.pauli):
        return Fraction(0)
    return term.coeff * (-1 if sum(bit for c, bit in zip(term.pauli, state) if c == 'Z') % 2 else 1)


def certify_hamiltonian(terms: Iterable[PauliTerm], *, state: Sequence[int] | None = None,
                        omitted_coefficients: Iterable[Fraction] = (), sqrt_digits=12):
    terms, state = _validated(terms, state)
    groups = []; comparisons = 0
    nonidentity = [t for t in terms if any(c != 'I' for c in t.pauli)]
    identity = sum((t.coeff for t in terms if all(c == 'I' for c in t.pauli)), Fraction(0))
    for term in nonidentity:
        for group in groups:
            fits = True
            for other in group:
                comparisons += 1
                if not _anti(term.pauli, other.pauli):
                    fits = False
                    break
            if fits:
                group.append(term)
                break
        else:
            groups.append([term])
    lower = identity; records = []
    for group in groups:
        squared = sum((t.coeff ** 2 for t in group), Fraction(0))
        lo, hi = _sqrt_interval(squared, sqrt_digits)
        lower -= hi
        records.append({'terms': [{'pauli': t.pauli, 'coeff': str(t.coeff)} for t in group],
                        'sum_coeff_squared': str(squared), 'sqrt_lower': str(lo), 'sqrt_upper': str(hi)})
    omitted = [Fraction(x) for x in omitted_coefficients]
    eta = sum(map(abs, omitted), Fraction(0))
    upper = sum((_expectation(t, state) for t in terms), Fraction(0))
    return {'lower': str(lower - eta), 'upper': str(upper + eta),
            'identity': str(identity), 'omitted_norm': str(eta),
            'omitted_coefficients': list(map(str, omitted)), 'state': list(state),
            'cliques': records,
            'search': {'method': 'greedy_anticommuting_partition', 'terms': len(nonidentity),
                       'cliques': len(groups), 'anticommutation_comparisons': comparisons,
                       'sqrt_decimal_digits': sqrt_digits}}


def check_certificate(terms, certificate):
    """Verify the supplied witness without rerunning partition or sqrt search.

The omitted norm is an assumption carried by the certificate. Its justification
must come from a separate parameter/error certificate; this checker cannot infer
unknown physics from an asserted error budget.
"""
    try:
        terms, state = _validated(terms, certificate['state'])
        width = len(terms[0].pauli)
        identity = sum((t.coeff for t in terms if all(c == 'I' for c in t.pauli)), Fraction(0))
        if Fraction(certificate['identity']) != identity:
            return False
        expected = Counter((t.pauli, t.coeff) for t in terms if any(c != 'I' for c in t.pauli))
        seen = Counter(); lower = identity
        for group in certificate['cliques']:
            entries = [PauliTerm(t['pauli'], Fraction(t['coeff'])) for t in group['terms']]
            if not entries or any(len(t.pauli) != width for t in entries):
                return False
            if any(not _anti(t.pauli, u.pauli) for i, t in enumerate(entries) for u in entries[i+1:]):
                return False
            squared = sum((t.coeff ** 2 for t in entries), Fraction(0))
            lo, hi = Fraction(group['sqrt_lower']), Fraction(group['sqrt_upper'])
            if Fraction(group['sum_coeff_squared']) != squared or not (0 <= lo <= hi and lo ** 2 <= squared <= hi ** 2):
                return False
            lower -= hi
            seen.update((t.pauli, t.coeff) for t in entries)
        if seen != expected:
            return False
        eta = sum((abs(Fraction(v)) for v in certificate['omitted_coefficients']), Fraction(0))
        if Fraction(certificate['omitted_norm']) != eta:
            return False
        upper = sum((_expectation(t, state) for t in terms), Fraction(0))
        return Fraction(certificate['lower']) == lower - eta and Fraction(certificate['upper']) == upper + eta
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError, AttributeError, OverflowError):
        return False
