"""Exact scalar bounds for one Hubbard bond.

The formula is for all local occupations (0 through 4 electrons), so it does
not assume a globally fixed two-electron bond sector.
"""
from fractions import Fraction
from math import isqrt
from research.constructive_response_20260916.composition_exact import is_psd


def sqrt_up(x: Fraction, scale: int = 10**12) -> Fraction:
    if x < 0 or not isinstance(scale,int) or scale <= 0:
        raise ValueError("nonnegative radicand and positive integer scale required")
    n = x.numerator * scale * scale // x.denominator
    r = Fraction(isqrt(n), scale)
    while r * r < x:
        r += Fraction(1, scale)
    return r


def bond_b_upper(a: Fraction, t: Fraction = Fraction(1)) -> Fraction:
    """Certified (possibly slightly conservative) b with T+aD+bI >= 0.

    For a >= 0 the exact optimum is max(t, (sqrt(a^2+16t^2)-a)/2).
    ``sqrt_up`` makes the returned rational safely large.
    """
    if a < 0 or t < 0:
        raise ValueError("this bound requires a,t >= 0")
    radical = sqrt_up(a*a + 16*t*t)
    return max(t, (radical-a)/2)


def bond_spectrum_lower(a: Fraction, t: Fraction = Fraction(1)) -> Fraction:
    """Exact symbolic lower eigenvalue represented by a certified radical bound."""
    return -bond_b_upper(a, t)


def graph_bound(U: Fraction, t: Fraction, z: int, a: Fraction) -> tuple[Fraction, Fraction]:
    """Return (effective doublon coefficient, extensive constant).

    With each edge assigned once, sum_i d_i appears at most z times, so
    T >= -a*z*D - b*|E|.  This is an intentionally transparent overlap bound.
    """
    if z < 0 or U < 0:
        raise ValueError("invalid graph parameters")
    b = bond_b_upper(a, t)
    return U-a*z, b


def _add(A, B):
    return [[A[i][j]+B[i][j] for j in range(len(A))] for i in range(len(A))]


def _hop(m, n):
    """Canonical c_n^dagger c_m on four ordered spin modes."""
    M = [[Fraction(0) for _ in range(16)] for _ in range(16)]
    for s in range(16):
        if (s >> n) & 1 or not ((s >> m) & 1): continue
        sign1 = -1 if (s & ((1 << m)-1)).bit_count() & 1 else 1
        s1 = s ^ (1 << m)
        sign2 = -1 if (s1 & ((1 << n)-1)).bit_count() & 1 else 1
        s2 = s1 | (1 << n)
        M[s2][s] += sign1*sign2
    return M


def bond_matrix(a, t=Fraction(1), b=None):
    """Independent CAR construction of T+aD+bI on ALL local occupations."""
    if a < 0 or t < 0: raise ValueError("this bound requires a,t >= 0")
    T = _add(_hop(0,2), _hop(2,0))
    T = _add(T, _add(_hop(1,3), _hop(3,1)))
    D = [[Fraction(0) for _ in range(16)] for _ in range(16)]
    for s in range(16):
        D[s][s] = int((s&3)==3) + int((s&12)==12)
    if b is None:
        b = bond_b_upper(a,t)
    return [[-t*T[i][j] + a*D[i][j] + (b if i == j else 0)
             for j in range(16)] for i in range(16)]


def verify_bond_all_sectors(a, t=Fraction(1)):
    b=bond_b_upper(a,t)
    return is_psd(bond_matrix(a,t,b)),b
