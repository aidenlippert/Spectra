"""Nonaccepting numerical adapter for two previously verified residual directions."""
from fractions import Fraction as F
from residual_coherence_obstruction import pure_coherence
LABELS = {'358,409,1':{358:1,409:1}, '103,358,1':{103:1,358:1}}
def actions(terms):
    if type(terms) is not dict or any(k not in LABELS for k in terms):
        raise ValueError('Canonical residual direction required')
    result=[{} for _ in range(4096)]
    for label,coefficient in terms.items():
        denominator,_,matrix,_=pure_coherence(LABELS[label])
        for (r,c),a in matrix.items():
            result[c][r]=result[c].get(r,F(0))+F(coefficient)*F(a,denominator)
    return result
