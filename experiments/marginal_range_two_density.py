"""Bounded exact next-nearest density profiles and diagonal actions."""
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact


def local_profile(sites, coupling, values=None):
    if type(sites) is not int or sites not in (4, 6):
        raise ValueError('Range-two local profile requires four or six sites')
    coupling = _exact(coupling)
    if values is None:
        values = [F(sites-1, sites-2)*coupling]*(sites-2)
    if type(values) is not list or len(values) != sites-2:
        raise ValueError('Range-two density profile has wrong length')
    values = [_exact(v) for v in values]
    if values != values[::-1] or sum(values) != (sites-1)*coupling:
        raise ValueError('Range-two profile must preserve reflection and translated mean')
    return values


def diagonal_value(state, weights):
    # q_i = n_{i,up}+n_{i,down}-1. Validated callers supply L-2 weights.
    charges = [((state >> (2*i)) & 3).bit_count()-1 for i in range(len(weights)+2)]
    return sum((w*charges[i]*charges[i+2] for i, w in enumerate(weights)), F(0))
