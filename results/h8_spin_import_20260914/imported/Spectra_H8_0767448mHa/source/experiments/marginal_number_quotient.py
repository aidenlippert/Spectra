"""Degree-preserving exact number-slice reduction with explicit ideal witnesses.

This is coefficient algebra, not positivity or state enumeration. A fixed
maximum degree controls the monomial basis; growing degree can still make it
combinatorial. The existing positivity checker remains the acceptance gate.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb


def _add(target, source, scale=F(1)):
    for mask, value in source.items():
        updated = target.get(mask, F(0))+scale*value
        if updated:
            target[mask] = updated
        else:
            target.pop(mask, None)


def _order(mask):
    return mask.bit_count(), mask


class NumberSliceQuotient:
    def __init__(self, variables, population, degree, monomial_budget=4096):
        if (type(variables) is not list or not variables
                or len(set(variables)) != len(variables)
                or any(type(i) is not int or not 0 <= i < 64 for i in variables)):
            raise ValueError('Distinct variable indices in0..63 required')
        m = len(variables)
        if (type(population) is not int or not 0 < population < m
                or type(degree) is not int or not 0 <= degree <= min(population, m-population)):
            raise ValueError('Degree must be within both slice populations')
        count = sum(comb(m, k) for k in range(degree+1))
        if type(monomial_budget) is not int or monomial_budget < 1 or count > monomial_budget:
            raise ValueError('Coefficient monomial budget exceeded')
        self.variables = sorted(variables)
        self.mask = sum(1 << i for i in variables)
        self.population, self.degree = population, degree
        self.monomials = [sum(1 << i for i in chosen)
                          for k in range(degree+1)
                          for chosen in combinations(self.variables, k)]
        self.rows = {}
        for support in self.monomials:
            if support.bit_count() >= degree:
                continue
            row = {support: F(support.bit_count()-population)}
            row.update({support | (1 << i): F(1) for i in self.variables if not support & (1 << i)})
            witness = {support: F(1)}
            while row:
                pivot = max(row, key=_order)
                coefficient = row[pivot]
                if pivot not in self.rows:
                    self.rows[pivot] = ({k: v/coefficient for k, v in row.items()},
                                        {k: v/coefficient for k, v in witness.items()})
                    break
                prior, prior_witness = self.rows[pivot]
                _add(row, prior, -coefficient)
                _add(witness, prior_witness, -coefficient)
            else:
                raise ValueError('Unexpected dependent number relation')
        self.pivots = sorted(self.rows, key=_order, reverse=True)
        self.cache = {}
        self.stats = {'variables': m, 'population': population, 'maximum_degree': degree,
                      'monomials': count, 'relations': len(self.rows),
                      'quotient_dimension': count-len(self.rows), 'configuration_evaluations': 0}
        if self.stats['quotient_dimension'] != comb(m, degree):
            raise ValueError('Unexpected bounded-slice quotient rank')

    def reduce(self, polynomial):
        if any(type(mask) is not int or mask < 0 or mask & ~self.mask
               or mask.bit_count() > self.degree for mask in polynomial):
            raise ValueError('Polynomial outside this coefficient space')
        current = {mask: F(value) for mask, value in polynomial.items() if value}
        multiplier = {}
        for pivot in self.pivots:
            coefficient = current.get(pivot)
            if coefficient:
                row, witness = self.rows[pivot]
                _add(current, row, -coefficient)
                _add(multiplier, witness, coefficient)
        return current, multiplier

    def monomial(self, mask):
        if mask not in self.cache:
            self.cache[mask] = self.reduce({mask: F(1)})
        # Callers must not mutate these cached coefficient dictionaries.
        return self.cache[mask]


def complete_bounded_number_ideals(ring, polynomial):
    """Return remainder,Lalpha,Lbeta with no increase in occupation degree.

    The identity is in the same population-pruned Boolean ring as JointPolynomial.
    No full-population monomial lifting, physical state list or action is used.
    """
    current = {mask: F(value) for mask, value in polynomial.items() if value}
    ideals = []
    receipts = []
    for spin in ring.spin_masks:
        variables = [i for i in range(ring.modes) if spin & (1 << i)]
        degree = max(((mask & spin).bit_count() for mask in current), default=0)
        quotient = NumberSliceQuotient(variables, ring.target, degree)
        reduced, multiplier = {}, {}
        for mask, coefficient in current.items():
            remainder, witness = quotient.monomial(mask & spin)
            other = mask & ~spin
            _add(reduced, {key | other: value for key, value in remainder.items()}, coefficient)
            _add(multiplier, {key | other: value for key, value in witness.items()}, coefficient)
        current = reduced
        ideals.append(multiplier)
        receipts.append(quotient.stats)
    return current, ideals, receipts
