"""Exact bounded-degree singlet functional by occupation counting and spin twirl."""
from fractions import Fraction as F
from math import comb
from experiments.marginal_symbolic import canonical
from research.sector_quotient_20260914.fast_twirl import twirl


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def singlet_dimension(modes, particles):
    if type(modes) is not int or type(particles) is not int or modes % 2 or particles % 2 or not 0 <= particles <= modes:
        raise ValueError('Paired spin orbitals and even particle number required')
    s, n = modes//2, particles//2
    return choose(s, n)**2-choose(s, n+1)*choose(s, n-1)


def magnetic_trace(poly, modes, alpha, beta):
    """Unnormalized trace; only identical creation/annihilation sets contribute."""
    s = modes//2
    total = F(0)
    for w, value in canonical(poly).items():
        if any(type(i) is not int or not 0 <= i < modes for c, i in w):
            raise ValueError('Operator mode outside the supplied model')
        creators = tuple(i for c, i in w if c)
        annihilators = tuple(i for c, i in w if not c)
        if creators != annihilators:
            continue
        degree = len(creators)
        a = sum(i % 2 == 0 for i in creators)
        b = degree-a
        total += value*(-1)**(degree*(degree-1)//2)*choose(s-a, alpha-a)*choose(s-b, beta-b)
    return total


def singlet_trace(poly, modes, particles):
    dimension = singlet_dimension(modes, particles)
    p = canonical(poly)
    if any(len(w) > 6 or len(w) % 2 or sum(2*c-1 for c, i in w) for w in p):
        raise ValueError('Balanced polynomial through degree six required')
    averaged = twirl(p)
    n = particles//2
    return (magnetic_trace(averaged, modes, n, n)-magnetic_trace(averaged, modes, n+1, n-1))/dimension


def spatial_charge(poly):
    signatures = set()
    for w in poly:
        counts = {}
        for c, i in w:
            counts[i//2] = counts.get(i//2, 0)+2*c-1
        signatures.add(tuple((i, v) for i, v in sorted(counts.items()) if v))
    if len(signatures) != 1:
        raise ValueError('A sparse trace block requires a single spatial charge signature')
    return next(iter(signatures))
