"""Exact product Hubbard-dimer diagnostics for global elimination.

This is a structural control: it proves how a locally gapped eliminated
sector can acquire exponentially small global overlap.  It is not a
connected molecular certificate.
"""
from fractions import Fraction as F
from math import isqrt
import sys

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

U = F(8)
t = F(1)
RAD = F(80)


def sqrt_interval(x, scale=1 << 80):
    if x < 0 or scale <= 0:
        raise ValueError("positive rational and scale required")
    n = (x.numerator * scale * scale) // x.denominator
    lo = isqrt(n)
    while (lo + 1) * (lo + 1) * x.denominator <= x.numerator * scale * scale:
        lo += 1
    while lo * lo * x.denominator > x.numerator * scale * scale:
        lo -= 1
    hi = lo + (lo * lo * x.denominator != x.numerator * scale * scale)
    return F(lo, scale), F(hi, scale)


def dimer_energy_interval():
    lo, hi = sqrt_interval(RAD)
    return F(4) - hi / 2, F(4) - lo / 2


def bare_weight_interval():
    lo, hi = sqrt_interval(RAD)
    return (F(1, 2) + F(8, 2) / hi, F(1, 2) + F(8, 2) / lo)


def dressed_parameters():
    r = F(116434, 1000000)
    return (1-r*r)/(1+r*r), 2*r/(1+r*r)


def dressed_weight_interval(c=None, s=None):
    if c is None and s is None:
        c, s = dressed_parameters()
    elif c is None or s is None:
        raise ValueError("provide both rotation entries")
    if c * c + s * s != 1:
        raise ValueError("dressing must be an exactly normalized rational rotation")
    lo, hi = sqrt_interval(RAD)
    # Ground-state density of [[0,-2],[-2,8]] in the {|S>,|D>} basis.
    num = U * (c * c - s * s) / 2 + 4 * c * s
    endpoints=[F(1, 2)+num/hi,F(1, 2)+num/lo]
    return min(endpoints),max(endpoints)


def powers(interval, L):
    if type(L) is not int or L < 1:
        raise ValueError("positive integer family size required")
    lo, hi = interval
    if not 0 <= lo <= hi <= 1:
        raise ValueError("probability enclosure in [0,1] required")
    return lo ** L, hi ** L


def diagnostics(L, dressed=False):
    weight = dressed_weight_interval() if dressed else bare_weight_interval()
    e_lo, e_hi = dimer_energy_interval()
    p_lo, p_hi = powers(weight, L)
    # E0=L*e. For the dressed projector the retained expectation is not zero.
    e0_lo, e0_hi = L * e_lo, L * e_hi
    if not 0 < p_lo <= p_hi < 1:
        raise ValueError("strictly interior overlap needed for this obstruction")
    if dressed:
        c, s = dressed_parameters()
        local_a = -4*c*s + 8*s*s
        gap_hi = L * (local_a - e_lo) * p_hi / (1 - p_hi)
    else:
        gap_hi = (-e0_lo) * p_hi / (1 - p_hi)
    response_lo = sqrt_interval((1 - p_hi) / p_hi)[0]
    # The full product Hamiltonian has a unique ground state and gap -E_d in
    # the balanced sector (one local m_s=0 triplet realizes that excitation).
    # Every normalized Q vector has ground-state overlap squared <=1-p.
    gap_lo=(-e_hi)*p_lo
    return {
        'L': L, 'dressed': dressed,
        'local_energy_interval_over_t': [str(e_lo), str(e_hi)],
        'global_energy_interval_over_t': [str(e0_lo), str(e0_hi)],
        'local_retained_weight_interval': [str(weight[0]), str(weight[1])],
        'global_overlap_interval': [str(p_lo), str(p_hi)],
        'global_Qgap_upper_over_t': str(gap_hi),
        'global_Qgap_lower_over_t': str(gap_lo),
        'response_norm_lower': str(response_lo),
        'global_states_constructed': 0,
    }


def run():
    rows = []
    for dressed in (False, True):
        rows.append({'family': 'rational_dressed' if dressed else 'bare',
                     'rows': [diagnostics(L, dressed) for L in
                              (1, 2, 4, 8, 16, 32, 64, 128, 256)]})
    return {'kind': 'product_hubbard_dimer_overlap_v1', 'U': str(U), 't': str(t),
            'local_hamiltonian': [["0", "-2"], ["-2", "8"]],
            'dressing_ratio': '116434/1000000',
            'families': rows, 'scope': 'Disconnected product control; no connected-proof or physical-scaling claim.'}


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2))
