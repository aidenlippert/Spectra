"""Compact reflection-odd quadratic charge operators on five sites."""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact


PAIRS = tuple((i, j) for i in range(5) for j in range(i, 5) if (i, j) < (4-j, 4-i))
LABELS = {f'{i},{j}': (i, j) for i, j in PAIRS}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= len(PAIRS):
        raise ValueError('One through six compact quadratic components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical reflection-odd quadratic pair required')
    result = {key: _exact(value) for key, value in source.items()}
    result = {key: value for key, value in result.items() if value}
    if not result:
        raise ValueError('Nonzero compact quadratic telescope required')
    return result


def five_site_value(state, terms):
    q = [((state >> (2*i)) & 3).bit_count()-1 for i in range(5)]
    return sum((a*(q[LABELS[key][0]]*q[LABELS[key][1]]-q[4-LABELS[key][1]]*q[4-LABELS[key][0]])
                for key, a in terms.items()), F(0))


def local_value(state, terms):
    return five_site_value(state & 1023, terms)-five_site_value(state >> 2, terms)
