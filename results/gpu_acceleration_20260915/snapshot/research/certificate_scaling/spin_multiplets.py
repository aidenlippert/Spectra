"""Exact rational SU(2) multiplets and positive invariant-square controls.

These are algebra primitives, not a completed spin-adapted discovery backend.
"""
from fractions import Fraction as F
from math import comb,lcm
from experiments.marginal_symbolic import canonical,validate_word,add,scale,product
from experiments.marginal_hunt_car import adj
from experiments.marginal_spin_reduction import spin_operators,commutator


def multiplet(highest,modes,two_spin):
    if type(two_spin) is not int or two_spin not in (0,1,2,3):
        raise ValueError('Supported spins are 0, 1/2, 1, 3/2')
    plus,minus,z=spin_operators(modes)
    if not isinstance(highest,dict) or not highest:
        raise ValueError('Nonzero highest-weight polynomial required')
    for w,c in highest.items():
        validate_word(w,modes,3)
        if type(c) not in (int,F):raise ValueError('Exact real rational coefficients required')
    p=canonical(highest)
    charges={sum(2*c-1 for c,_ in w) for w in p}
    if not p or len(charges)!=1:raise ValueError('One nonzero number-charge required')
    if commutator(plus,p) or commutator(z,p)!=scale(p,F(two_spin,2)):
        raise ValueError('Highest-weight identities failed')
    descendants=[p]
    for r in range(1,two_spin+1):
        descendants.append(scale(commutator(minus,descendants[-1]),F(1,r)))
    if commutator(minus,descendants[-1]):raise ValueError('Multiplet did not terminate')
    for r,q in enumerate(descendants):
        if commutator(z,q)!=scale(q,F(two_spin,2)-r):raise ValueError('Weight identity failed')
        if r and commutator(plus,q)!=scale(descendants[r-1],two_spin-r+1):
            raise ValueError('Raising identity failed')
    denominator=lcm(*(comb(two_spin,r) for r in range(two_spin+1)))
    multiplicities=[denominator//comb(two_spin,r) for r in range(two_spin+1)]
    return descendants,multiplicities


def invariant_square(highest,modes,two_spin):
    descendants,weights=multiplet(highest,modes,two_spin)
    return add(*(scale(product(canonical(adj(p)),p),weight) for p,weight in zip(descendants,weights)))


def export_invariant_square(highest,modes,two_spin):
    descendants,weights=multiplet(highest,modes,two_spin)
    denominator=lcm(*(c.denominator for p in descendants for c in p.values()))
    blocks=[]
    for p,weight in zip(descendants,weights):
        words=sorted(p,key=lambda w:(len(w),w))
        row=[int(p[w]*denominator) for w in words]
        blocks.append({'words':words,'factor':[row[:] for _ in range(weight)]})
    return blocks,denominator
