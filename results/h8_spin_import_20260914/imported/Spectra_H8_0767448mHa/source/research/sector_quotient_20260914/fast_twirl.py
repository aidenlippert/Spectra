"""Exact proposal-side SU(2) projection on already canonical CAR words.

Spin partners are adjacent indices. Flipping one index cannot change the order
within a sorted creation/annihilation group: it either keeps that order or
duplicates an index and gives zero. This avoids repeated general CAR reduction.
The accepting molecular replay continues to use the sealed original projector.
"""
from fractions import Fraction as F
from experiments.marginal_symbolic import add, scale
from research.certificate_scaling.spin_basis import spin_weight2


def ladder(poly, raising):
    out = {}
    for w, value in poly.items():
        present = set(w)
        for k, (creation, i) in enumerate(w):
            if not ((creation == 1 and i % 2 == (1 if raising else 0)) or
                    (creation == 0 and i % 2 == (0 if raising else 1))):
                continue
            replacement = (creation, i ^ 1)
            if replacement in present: continue
            v = w[:k]+(replacement,)+w[k+1:]
            out[v] = out.get(v, F(0))+(value if creation else -value)
    return {w: value for w, value in out.items() if value}


def casimir(poly):
    z = {w: value*F(spin_weight2(w)**2, 4) for w, value in poly.items() if spin_weight2(w)}
    return add(z, scale(add(ladder(ladder(poly, False), True), ladder(ladder(poly, True), False)), F(1, 2)))


def twirl(poly):
    for w, value in poly.items():
        left = tuple((c, i) for c, i in w if c)
        right = tuple((c, i) for c, i in w if not c)
        if (len(w) > 6 or len(w) % 2 or w != left+right or
            left != tuple(sorted(set(left))) or right != tuple(sorted(set(right))) or
            any(c not in (0, 1) or type(i) is not int or i < 0 for c, i in w)):
            raise ValueError('Expected an even canonical CAR polynomial through degree six')
        if type(value) not in (int, F): raise ValueError('Exact real rational coefficients required')
    out = dict(poly)
    for j in range(1, max((len(w)//2 for w in out), default=0)+1):
        out = add(out, scale(casimir(out), F(-1, j*(j+1))))
    if casimir(out): raise AssertionError('Exact spin projection failed Casimir check')
    return out
