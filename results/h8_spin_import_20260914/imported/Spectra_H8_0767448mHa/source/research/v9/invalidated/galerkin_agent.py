"""Bounded trajectory-wide residual driven adaptive Galerkin baseline.

This is a conventional exact arithmetic proposer.  The checker remains v7:
the projected recurrence is exported as an ordinary polynomial and every
residual (including omitted terms) is reconstructed independently.
"""
from collections import defaultdict
from fractions import Fraction as F

from .v7_certificate import (BudgetExceeded, Piece, clean, derive_certificate,
                             rational, check_certificate, norm_witness)


def projected_taylor(gen, initial, basis, duration, max_order=24):
    """Return the exact projected Taylor coefficients and omitted residuals.

    ``basis`` is fixed for this call.  Residuals are (I-P)G c_k, with every
    coefficient scanned, so collisions and signs are retained.
    """
    if type(max_order) is not int or not 0 <= max_order <= 24:
        raise ValueError('order cap')
    t = rational(duration)
    if t <= 0: raise ValueError('positive duration required')
    b = set(basis)
    if not b: raise ValueError('nonempty basis')
    if any(not isinstance(p, str) or len(p) != gen.n or any(x not in 'IXYZ' for x in p) for p in b):
        raise ValueError('invalid basis label')
    current = {p: c for p, c in clean(initial, gen.n, gen.max_terms).items() if p in b}
    coeffs = [current]
    omitted = []
    for k in range(max_order + 1):
        full = gen.apply(current)
        omitted.append({p: v for p, v in full.items() if p not in b and v})
        if k == max_order:
            break
        nxt = {p: v / F(k + 1) for p, v in full.items() if p in b and v}
        coeffs.append(clean(nxt, gen.n, gen.max_terms))
        current = nxt
    return tuple(coeffs), tuple(omitted)


def _score(omitted, duration):
    t = rational(duration)
    score = defaultdict(F)
    for k, op in enumerate(omitted):
        weight = t ** (k + 1) / F(k + 1)
        for p, value in op.items():
            score[p] += abs(value) * weight
    return dict(score)


def adaptive_galerkin(gen, initial, duration, tolerance, *, max_order=24,
                      max_expansions=16, batch=None, grouping='l1',
                      integration_basis='power'):
    """Construct a bounded exact Galerkin candidate or an honest refusal.

    Basis growth is deterministic: labels are sorted by accumulated weighted
    omitted residual (then Pauli label).  Every restart and failed round is
    retained in the returned work report.
    """
    tol = rational(tolerance); t = rational(duration)
    if tol < 0 or t <= 0: raise ValueError('duration/tolerance')
    seed = clean(initial, gen.n, gen.max_terms)
    basis = set(seed)
    if not basis: raise ValueError('nonempty initial support required')
    if type(max_expansions) is not int or not 0 <= max_expansions <= 16:
        raise ValueError('expansion cap')
    if batch is None: batch = max(4, len(basis)//2)
    if type(batch) is not int or not 1 <= batch <= 512: raise ValueError('batch cap')
    rounds = []
    for expansion in range(max_expansions + 1):
        before = len(gen.cache)
        # Increment order until the in-basis final row is small enough.  This
        # keeps Taylor tail and omitted residual decisions separate.
        chosen = max_order
        for order in range(max_order + 1):
            probe, probe_omitted = projected_taylor(gen, seed, basis, t, order)
            inside = gen.apply(probe[-1])
            tail = sum((abs(v) * t ** (order + 1) / F(order + 1)
                        for p, v in inside.items() if p in basis), F(0))
            if tail <= tol / 4:
                chosen = order
                break
        coeffs, omitted = projected_taylor(gen, seed, basis, t, chosen)
        piece = Piece(t, coeffs)
        cert = derive_certificate(gen, seed, [piece], grouping, integration_basis)
        independent = type(gen)(dict(gen.h), gen.gamma, gen.n, gen.max_terms)
        replay = check_certificate(independent, seed, [piece], cert, tol, expected_time=t)
        scores = _score(omitted, t)
        candidates = sorted((p for p in scores if p not in basis),
                            key=lambda p: (-scores[p], p))
        tail = sum((abs(v) * t ** (chosen + 1) / F(chosen + 1)
                    for p, v in (gen.apply(coeffs[-1]) if coeffs else {}).items()
                    if p in basis), F(0))
        rounds.append({'round': expansion, 'basis_size': len(basis),
                       'new_columns': len(gen.cache) - before,
                       'omitted_score': {p: str(scores[p]) for p in candidates},
                       'tail_bound': str(tail), 'status': replay['status'],
                       'bound': replay.get('bound')})
        if replay['status'] == 'certified':
            return {'status': 'certified', 'piece': piece, 'certificate': cert,
                    'basis': tuple(sorted(basis)), 'rounds': rounds,
                    'expansions': expansion, 'checker': replay,
                    'generator_cost': dict(gen.cost)}
        if expansion >= max_expansions or not candidates:
            break
        add = candidates[:min(batch, 512 - len(basis))]
        basis.update(add)
        if len(basis) > gen.max_terms: raise BudgetExceeded('basis term budget')
    return {'status': 'refused', 'reason': 'bounded Galerkin expansion exhausted',
            'basis': tuple(sorted(basis)), 'rounds': rounds,
            'expansions': len(rounds) - 1, 'generator_cost': dict(gen.cost)}


# Short aliases used by experiments and downstream tests.
galerkin = adaptive_galerkin
