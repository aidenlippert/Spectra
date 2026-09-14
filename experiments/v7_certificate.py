"""Exact, model-conditional finite-time observable residual certificates.

A candidate is a piecewise polynomial with rational Pauli coefficients. Its
residual is recomputed; an untrusted proposer supplies only norm partitions.
Uniform local depolarization is a declared physical assumption, not a tunable
certificate parameter. No physical applicability claim follows from acceptance.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
import re
from math import comb
from .pauli import commutator_i
from .certificates import _anti, _sqrt_interval


class BudgetExceeded(ValueError):
    pass


def rational(x):
    if not isinstance(x, (int, F)) or isinstance(x, bool):
        raise ValueError('exact rational coefficients required')
    x = F(x)
    if max(x.numerator.bit_length(), x.denominator.bit_length()) > 8192:
        raise BudgetExceeded('rational bit-length budget')
    return x


def proof_fraction(x):
    """Parse bounded exact witness text before allocating large integers."""
    if isinstance(x, str):
        if len(x) > 5002 or not re.fullmatch(r"-?[0-9]{1,2500}(?:/[0-9]{1,2500})?", x):
            raise ValueError("invalid or oversized rational witness")
        x = F(x)
    return rational(x)


def clean(op, n, cap=2048):
    if not isinstance(op, dict) or len(op) > cap:
        raise BudgetExceeded('Pauli term budget or invalid mapping')
    out = {}
    for p, c in op.items():
        if not isinstance(p, str) or len(p) != n or any(x not in 'IXYZ' for x in p):
            raise ValueError('invalid Pauli label')
        c = rational(c)
        if c: out[p] = c
    return out


def add(a, b, scale=F(1)):
    out = dict(a)
    for p, c in b.items():
        out[p] = out.get(p, F(0)) + scale*c
        if not out[p]: del out[p]
    return out


def evaluate(coefficients, t):
    out = {}
    for c in reversed(coefficients):
        out = add({p: v*t for p,v in out.items()}, c)
    return out


@dataclass(frozen=True)
class Piece:
    duration: F
    coefficients: tuple


class Generator:
    """Column action i[H,P] - gamma * weight(P) P with explicit work counts."""
    def __init__(self, h, gamma, n, max_terms=2048):
        if type(n) is not int or not 1 <= n <= 32:
            raise ValueError('1 <= n <= 32 required')
        if type(max_terms) is not int or not 1 <= max_terms <= 4096:
            raise ValueError('invalid term cap')
        self.n, self.max_terms = n, max_terms
        self.h = clean(h, n, 256)
        self.gamma = rational(gamma)
        if self.gamma < 0: raise ValueError('negative depolarization')
        self.cache = {}
        self.cost = dict(hamiltonian_pairs=0, coefficient_multiply_adds=0,
                         column_requests=0, cache_hits=0, peak_terms=0,
                         max_rational_bits=0, basis_multiply_adds=0)

    def apply(self, op):
        op = clean(op, self.n, self.max_terms)
        out = {}
        for p, c in op.items():
            self.cost['column_requests'] += 1
            if p not in self.cache:
                if len(self.cache) >= self.max_terms:
                    raise BudgetExceeded('cached column budget')
                self.cost['hamiltonian_pairs'] += len(self.h)
                col = commutator_i(self.h, p)
                rate = -self.gamma * sum(x != 'I' for x in p)
                if rate: col = add(col, {p: rate})
                self.cache[p] = col
            else: self.cost['cache_hits'] += 1
            for q, d in self.cache[p].items():
                self.cost['coefficient_multiply_adds'] += 1
                out[q] = out.get(q, F(0)) + c*d
                if not out[q]: del out[q]
            if len(out) > self.max_terms: raise BudgetExceeded('output term budget')
        out = clean(out, self.n, self.max_terms)
        self.cost['peak_terms'] = max(self.cost['peak_terms'], len(op), len(out))
        for c in out.values():
            self.cost['max_rational_bits'] = max(self.cost['max_rational_bits'], c.numerator.bit_length(), c.denominator.bit_length())
        return out


def residual_records(gen, initial, pieces, integration_basis='power'):
    if integration_basis not in ('power', 'bernstein'):
        raise ValueError('invalid integration basis')
    """Return exact (polynomial/jump label, operator, integration weight)."""
    initial = clean(initial, gen.n, gen.max_terms)
    if not isinstance(pieces, (list, tuple)) or not 1 <= len(pieces) <= 32:
        raise ValueError('1..32 pieces required')
    count = 0
    normalized = []
    for piece in pieces:
        if not isinstance(piece, Piece): raise ValueError('Piece required')
        t = rational(piece.duration)
        if t <= 0: raise ValueError('positive duration required')
        cs = piece.coefficients
        if not isinstance(cs, (list, tuple)) or not 1 <= len(cs) <= 33:
            raise ValueError('polynomial degree <=32 required')
        cs = tuple(clean(c, gen.n, gen.max_terms) for c in cs)
        count += sum(len(c) for c in cs)
        if count > 50000: raise BudgetExceeded('total coefficient budget')
        normalized.append((t, cs))
    records = []
    previous = initial
    for j, (t, cs) in enumerate(normalized):
        records.append((f'jump:{j}', add(cs[0], previous, F(-1)), F(1)))
        residuals = []
        for k, c in enumerate(cs):
            derivative = {p:(k+1)*v for p,v in cs[k+1].items()} if k+1 < len(cs) else {}
            residuals.append(add(derivative, gen.apply(c), F(-1)))
        while len(residuals)>1 and not residuals[-1]: residuals.pop()
        m = len(residuals)-1
        for k in range(m+1):
            if integration_basis == 'power':
                op, weight = residuals[k], t**(k+1)/F(k+1)
            else:
                op = {}
                for ell in range(k+1):
                    factor = t**ell * F(comb(k,ell), comb(m,ell))
                    gen.cost['basis_multiply_adds'] += len(residuals[ell])
                    op = add(op, residuals[ell], factor)
                op, weight = clean(op,gen.n,gen.max_terms), t/F(m+1)
            records.append((f'residual:{j}:{k}', op, weight))
        previous = clean(evaluate(cs, t), gen.n, gen.max_terms)
    return records


def norm_witness(op, grouping='l1'):
    """Conventional partitions, with actual pair-comparison counts."""
    if grouping not in ('l1', 'firstfit', 'weighted'):
        raise ValueError('unsupported grouping')
    order = sorted(op, key=(lambda p:(-abs(op[p]), p)) if grouping == 'weighted' else None)
    groups, comparisons = [], 0
    for p in order:
        accepted = False
        if grouping != 'l1':
            for group in groups:
                fits = True
                for q in group:
                    comparisons += 1
                    if not _anti(p, q): fits = False; break
                if fits:
                    group.append(p); accepted = True; break
        if not accepted: groups.append([p])
    result = []
    for group in groups:
        square = sum((op[p]**2 for p in group), F(0))
        hi = abs(op[group[0]]) if len(group) == 1 else _sqrt_interval(square, 16)[1]
        result.append(dict(labels=group, upper=str(hi)))
    return result, dict(group_comparisons=comparisons, sqrt_enclosures=sum(len(g)>1 for g in groups), sorted_terms=len(order))


def derive_certificate(gen, initial, pieces, grouping='l1', integration_basis='power'):
    records = residual_records(gen, initial, pieces, integration_basis)
    witnesses, construction = {}, dict(group_comparisons=0, sqrt_enclosures=0, sorted_terms=0)
    claimed = F(0)
    for label, op, weight in records:
        groups, cost = norm_witness(op, grouping)
        witnesses[label] = groups
        for k, v in cost.items(): construction[k] += v
        claimed += weight * sum((F(g['upper']) for g in groups), F(0))
    return dict(schema='v7-residual-1', integration_basis=integration_basis, witnesses=witnesses, claimed_bound=str(claimed),
                construction_cost=construction, generator_cost=dict(gen.cost))


def check_certificate(gen, initial, pieces, witness, tolerance, *, expected_time):
    """Recompute residuals and verify norm witnesses independently of search.

    Acceptance bounds every density-state expectation at the final time, and
    at intermediate times by the corresponding accumulated bound. A supplied
    negative/malformed certificate fails closed. Work counters survive failure.
    """
    cost = dict(pair_checks=0, squares=0, scalar_bound_ops=0)
    try:
        tolerance = rational(tolerance)
        if tolerance < 0: raise ValueError('negative tolerance')
        requested_time = rational(expected_time)
        if requested_time <= 0: raise ValueError('positive requested horizon required')
        records = residual_records(gen, initial, pieces, witness['integration_basis'])
        actual_time = sum((rational(p.duration) for p in pieces), F(0))
        if actual_time != requested_time: raise ValueError('requested horizon mismatch')
        if witness['schema'] != 'v7-residual-1': raise ValueError('schema')
        ws = witness['witnesses']
        if set(ws) != {label for label, _, _ in records}: raise ValueError('record coverage')
        bound = F(0)
        for label, op, weight in records:
            groups = ws[label]
            if not isinstance(groups, list) or len(groups) > len(op): raise ValueError('group count')
            seen = Counter(); norm = F(0)
            for group in groups:
                labels = group['labels']
                if not isinstance(labels, list) or not labels or len(labels) > len(op): raise ValueError('invalid group')
                for i, p in enumerate(labels):
                    if p not in op: raise ValueError('extra term')
                    for q in labels[i+1:]:
                        cost['pair_checks'] += 1
                        if not _anti(p,q): raise ValueError('non-anticommuting group')
                seen.update(labels)
                hi = proof_fraction(group['upper'])
                square = sum((op[p]**2 for p in labels), F(0))
                cost['squares'] += len(labels)+1
                if hi < 0 or hi*hi < square: raise ValueError('invalid norm upper bound')
                norm += hi
            if seen != Counter(op.keys()): raise ValueError('term coverage')
            bound += norm*weight
            cost['scalar_bound_ops'] += 2
        if bound != proof_fraction(witness['claimed_bound']): raise ValueError('claimed bound mismatch')
        return dict(status='certified' if bound <= tolerance else 'over_tolerance', bound=str(bound), time=str(actual_time),
                    checking_cost=cost, generator_cost=dict(gen.cost))
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, ZeroDivisionError, OverflowError) as exc:
        return dict(status='rejected', reason=str(exc), checking_cost=cost, generator_cost=dict(gen.cost))
