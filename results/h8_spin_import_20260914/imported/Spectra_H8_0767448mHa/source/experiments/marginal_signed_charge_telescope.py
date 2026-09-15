"""Disjoint pattern basis for PH-even reflection-odd five-site charge functions."""
from fractions import Fraction as F
from itertools import product
from experiments.marginal_local_hubbard_block import _exact


def negative(q):
    return tuple(-v for v in q)


PATTERNS = {','.join(map(str,q)):q for q in product((-1,0,1),repeat=5)
            if q==min(q,negative(q)) and q<min(q[::-1],negative(q[::-1]))}
ORBIT = {pattern:(key,sign) for key,q in PATTERNS.items()
         for sign,patterns in [(1,(q,negative(q))),(-1,(q[::-1],negative(q[::-1])))]
         for pattern in patterns}


def coefficients(source):
    if type(source) is not dict or not 1<=len(source)<=len(PATTERNS):
        raise ValueError('One through 52 canonical signed-charge patterns required')
    if any(type(key) is not str or key not in PATTERNS for key in source):
        raise ValueError('Canonical signed-charge pattern required')
    terms={key:_exact(value) for key,value in source.items()}
    terms={key:value for key,value in terms.items() if value}
    if not terms:raise ValueError('Nonzero signed-charge telescope required')
    return terms


def component(state):
    q=tuple(((state>>(2*i))&3).bit_count()-1 for i in range(5))
    return ORBIT.get(q)


def five_site_value(state,terms):
    entry=component(state)
    return F(0) if entry is None else entry[1]*terms.get(entry[0],F(0))


def local_value(state,terms):
    return five_site_value(state&1023,terms)-five_site_value(state>>2,terms)
