"""Compact reflection-odd pair correlations of empty/double-site indicators."""
from fractions import Fraction as F
from experiments.marginal_quadratic_charge_telescope import LABELS as QUADRATIC_LABELS, coefficients as quadratic_coefficients

LABELS = {key: pair for key, pair in QUADRATIC_LABELS.items() if pair[0] != pair[1]}


def coefficients(source):
    terms = quadratic_coefficients(source)
    if any(key not in LABELS for key in source):
        raise ValueError('Distinct-site charge-square pair required')
    return terms


def five_site_value(state, terms):
    square = [int(((state >> (2*i)) & 3) in (0, 3)) for i in range(5)]
    return sum((a*(square[LABELS[key][0]]*square[LABELS[key][1]]
                   -square[4-LABELS[key][1]]*square[4-LABELS[key][0]])
                for key, a in terms.items()), F(0))


def local_value(state, terms):
    return five_site_value(state & 1023, terms)-five_site_value(state >> 2, terms)
