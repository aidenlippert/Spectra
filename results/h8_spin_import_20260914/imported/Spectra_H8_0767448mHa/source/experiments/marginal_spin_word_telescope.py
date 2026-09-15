"""Complete120-dimensional PH/spin-even, reflection-odd diagonal five-site space."""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact, _reflection


def spin_flip(state):
    return sum(1 << (m ^ 1) for m in range(10) if state >> m & 1)


def _basis():
    remaining=set(range(1024));patterns={};lookup={}
    while remaining:
        seed=min(remaining)
        positive={seed,1023^seed,spin_flip(seed),1023^spin_flip(seed)}
        negative={_reflection(s,5)[0] for s in positive}
        orbit=positive|negative
        if not orbit<=remaining:raise ValueError('Spin-word orbits do not partition Fock words')
        remaining.difference_update(orbit)
        if positive&negative:continue
        label=str(seed);pattern={s:1 for s in positive}|{s:-1 for s in negative}
        for s,value in pattern.items():
            if pattern.get(1023^s)!=value or pattern.get(spin_flip(s))!=value or pattern.get(_reflection(s,5)[0])!=-value:
                raise ValueError('Spin-word orbit does not have the required symmetry')
            lookup[s]=(label,value)
        patterns[label]=pattern
    if len(patterns)!=120:raise ValueError('Incomplete diagonal spin-word basis')
    return patterns,lookup


PATTERNS,ORBIT=_basis()


def coefficients(source):
    if type(source) is not dict or not 1<=len(source)<=120:
        raise ValueError('One through120 canonical spin-word components required')
    if any(type(key) is not str or key not in PATTERNS for key in source):
        raise ValueError('Canonical spin-word label required')
    terms={key:_exact(value) for key,value in source.items()}
    terms={key:value for key,value in terms.items() if value}
    if not terms:raise ValueError('Nonzero exact spin-word telescope required')
    return terms


def component(state):
    return ORBIT.get(state)


def five_site_value(state,terms):
    entry=component(state)
    return F(0) if entry is None else entry[1]*terms.get(entry[0],F(0))


def local_value(state,terms):
    return five_site_value(state&1023,terms)-five_site_value(state>>2,terms)
