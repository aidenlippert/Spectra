"""Bounded rational evaluation of finite complex exponential-polynomial sums.

All disks use the l1 norm: ``|(x,y)| = |x| + |y|``.  Thus a result is a
midpoint pair and a rational radius which bounds the l1 distance to the exact
value.  No floating point values are used in the certificate.
"""
from __future__ import annotations

from fractions import Fraction
from math import gcd

MAX_BITS = 8192
MAX_TERMS = 50_000
MAX_DEGREE = 32
MAX_ITER = 100_000

Pair = tuple[Fraction, Fraction]


def _q(x) -> Fraction:
    if isinstance(x, bool) or not isinstance(x, (Fraction, int)):
        raise TypeError("inputs must be exact Fraction or int")
    return Fraction(x)


def _bits(x: Fraction) -> int:
    return max(x.numerator.bit_length(), x.denominator.bit_length())


def _check(x: Fraction) -> Fraction:
    if _bits(x) > MAX_BITS:
        raise ValueError("rational bit cap exceeded")
    return x


def _padd(a: Pair, b: Pair) -> Pair:
    return (_check(a[0] + b[0]), _check(a[1] + b[1]))


def _pscale(a: Pair, c: Fraction) -> Pair:
    return (_check(a[0] * c), _check(a[1] * c))


def _pnorm(a: Pair) -> Fraction:
    return abs(a[0]) + abs(a[1])


def _mul_disk(a: Pair, ar: Fraction, b: Pair, br: Fraction) -> tuple[Pair, Fraction]:
    # l1 is submultiplicative for complex multiplication.
    m = (_check(a[0] * b[0] - a[1] * b[1]),
         _check(a[0] * b[1] + a[1] * b[0]))
    r = _check(_pnorm(a) * br + _pnorm(b) * ar + ar * br)
    return m, r


def _exp_disk(z: Pair, tol: Fraction, cost: dict) -> tuple[Pair, Fraction]:
    """Compute exp(z) with an l1 disk, using scaling and exact Taylor terms."""
    mag = _pnorm(z)
    if mag == 0:
        return (Fraction(1), Fraction(0)), Fraction(0)
    # z/s <= 1/2 in l1, then square s times.
    s = 1
    while mag / s > Fraction(1, 2):
        s *= 2
        if s > MAX_ITER:
            raise ValueError("exponential scaling limit exceeded")
    w = _pscale(z, Fraction(1, s))
    # The omitted tail is bounded by the next term times 1/(1-r), where
    # r=|w|/(n+1) bounds all subsequent term ratios.
    term = (Fraction(1), Fraction(0))
    total = term
    n = 0
    rem = Fraction(3, 2)
    while True:
        n += 1
        term = _pscale((
            _check(term[0] * w[0] - term[1] * w[1]),
            _check(term[0] * w[1] + term[1] * w[0])), Fraction(1, n))
        total = _padd(total, term)
        # ratio of subsequent absolute Taylor terms is <= |w|/(n+1) <= 1/2.
        ratio = mag / (s * (n + 1))
        rem = _check(_pnorm(term) * ratio / (1 - ratio))
        cost["taylor_terms"] += 1
        if rem <= tol:
            break
        if n >= MAX_ITER:
            raise ValueError("Taylor iteration limit exceeded")
    radius = rem
    for _ in range(s.bit_length() - 1):
        total, radius = _mul_disk(total, radius, total, radius)
        cost["squarings"] += 1
        if radius > tol * 4 and s > 1 and _bits(radius) > MAX_BITS:
            raise ValueError("exponential certificate exceeded limits")
    return total, radius


def evaluate_mode_map(mapping, time, tolerance):
    """Evaluate several named sums while sharing exponent certificates."""
    T, tol = _q(time), _q(tolerance)
    if T < 0 or tol <= 0 or len(mapping) > 512:
        raise ValueError("invalid time, tolerance, or word count")
    words = {str(k): list(v) for k, v in mapping.items()}
    allrows = [row for rows in words.values() for row in rows]
    if len(allrows) > MAX_TERMS:
        raise ValueError("term count limit exceeded")
    groups = {}
    parsed = {}
    for name, rows in words.items():
        parsed[name] = []
        for row in rows:
            if len(row) != 5: raise ValueError("each term has five fields")
            a,b,d,cr,ci = row
            a,b,cr,ci = map(_q,(a,b,cr,ci))
            if a > 0 or type(d) is not int or d < 0 or d > MAX_DEGREE: raise ValueError("invalid term")
            item=(a,b,d,cr,ci); parsed[name].append(item)
            groups.setdefault((a,b), []).append((d,cr,ci,T))
    cost={"input_terms":len(allrows),"unique_exponents":len(groups),"taylor_terms":0,"squarings":0,"cache_hits":0}
    cache={}
    for key, rs in groups.items():
        weight=sum((abs(cr)+abs(ci))*T**d for d,cr,ci,_ in rs)
        local=tol/(max(1,len(words))*max(Fraction(1),weight)*16)
        while True:
            cache[key]=_exp_disk((_check(key[0]*T),_check(key[1]*T)),local,cost)
            if cache[key][1]*weight <= tol/max(1,len(words)): break
            local/=2
    result={}; total=Fraction(0)
    for name, rows in parsed.items():
        mid=(Fraction(0),Fraction(0)); er=Fraction(0)
        for a,b,d,cr,ci in rows:
            em,ee=cache[(a,b)]; coeff=(_check(cr*T**d),_check(ci*T**d))
            p,r=_mul_disk(coeff,Fraction(0),em,ee); mid=_padd(mid,p); er=_check(er+r)
        result[name]=mid; total=_check(total+er)
    return result,total,cost


def evaluate_modes(terms, time, tolerance):
    """Evaluate ``sum (cr+i ci) T**degree exp((a+i b)T)``.

    Returns ``((real, imag), absolute_error_bound, cost)``.  The bound is in
    l1 norm and includes every exponential and arithmetic propagation error.
    Inputs are required to be finite Fractions within documented caps.
    """
    T, tol = _q(time), _q(tolerance)
    if T < 0 or tol <= 0:
        raise ValueError("time must be nonnegative and tolerance positive")
    raw = list(terms)
    if len(raw) > MAX_TERMS:
        raise ValueError("term count limit exceeded")
    groups: dict[tuple[Fraction, Fraction], list[tuple[int, Fraction, Fraction]]] = {}
    for row in raw:
        if len(row) != 5:
            raise ValueError("each term has five fields")
        a, b, d, cr, ci = row
        a, b, cr, ci = map(_q, (a, b, cr, ci))
        if a > 0 or not isinstance(d, int) or d < 0 or d > MAX_DEGREE:
            raise ValueError("invalid exponent or polynomial degree")
        for x in (a, b, cr, ci, T, tol):
            _check(x)
        groups.setdefault((a, b), []).append((d, cr, ci))
    cost = {"input_terms": len(raw), "unique_exponents": len(groups),
            "taylor_terms": 0, "squarings": 0, "cache_hits": 0}
    out, err = (Fraction(0), Fraction(0)), Fraction(0)
    cache = {}
    for key, rows in groups.items():
        if key in cache:
            em, er = cache[key]; cost["cache_hits"] += 1
        else:
            weight = sum((abs(cr) + abs(ci)) * (T ** d) for d, cr, ci in rows)
            # Each group contributes at most weight*er; split the requested
            # budget over groups (and use one when the polynomial is zero).
            local_tol = tol / (max(1, len(groups)) * max(Fraction(1), weight) * 16)
            while True:
                em, er = _exp_disk((_check(key[0] * T), _check(key[1] * T)), local_tol, cost)
                if er * weight <= tol / max(1, len(groups)):
                    break
                local_tol /= 2
                if local_tol == 0:
                    raise ValueError("unable to meet requested tolerance")
            cache[key] = (em, er)
        for d, cr, ci in rows:
            tp = T ** d
            coeff = (_check(cr * tp), _check(ci * tp))
            mid, r = _mul_disk(coeff, Fraction(0), em, er)
            out = _padd(out, mid)
            err = _check(err + r)
    return out, err, cost
