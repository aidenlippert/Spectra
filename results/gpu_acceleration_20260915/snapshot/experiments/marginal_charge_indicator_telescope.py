"""Complete reflection-odd functions of five binary empty/double indicators."""
from fractions import Fraction as F
from itertools import combinations
from math import prod
from experiments.marginal_local_hubbard_block import _exact

ALL_LABELS = {','.join(map(str, sites)): sites for degree in range(1,5)
              for sites in combinations(range(5),degree)
              if sites < tuple(sorted(4-i for i in sites))}
LABELS = {key:sites for key,sites in ALL_LABELS.items() if len(sites)>=3}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= len(LABELS):
        raise ValueError('One through six higher indicator components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical triple or quadruple indicator subset required')
    terms = {key:_exact(value) for key,value in source.items()}
    terms = {key:value for key,value in terms.items() if value}
    if not terms:
        raise ValueError('Nonzero higher indicator telescope required')
    return terms


def five_site_value(state, terms):
    indicators = [int(((state >> (2*i)) & 3) in (0,3)) for i in range(5)]
    return sum((value*(prod(indicators[i] for i in ALL_LABELS[key])
                      -prod(indicators[4-i] for i in ALL_LABELS[key]))
                for key,value in terms.items()),F(0))


def local_value(state, terms):
    return five_site_value(state & 1023, terms)-five_site_value(state >> 2, terms)
