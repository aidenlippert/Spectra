"""Conventional exact rejection guards before Taylor norm proposals.

Copied from frozen V8 to preserve its evidence. The original checker and V8
source remain unchanged. This is supplied mathematics, not an acquired method.
"""
from fractions import Fraction as F
from math import gcd, factorial
from .v7_certificate import Piece, norm_witness, clean, rational, BudgetExceeded
from .pauli import commutator_i


def _lcm(a, b):
    return abs(a // gcd(a, b) * b)


def _denominator(*ops):
    d = 1
    for op in ops:
        for value in op.values():
            d = _lcm(d, F(value).denominator)
    return d


def guarded_taylor(gen, initial, duration, tolerance, max_order=24, guard="frobenius"):
    if guard not in ("max", "frobenius", "cascade"): raise ValueError("unknown norm guard")
    """Generate the same exact Taylor piece using integer recurrence numerators."""
    if type(max_order) is not int or not 0 <= max_order <= 24: raise ValueError('order cap')
    duration, tolerance = rational(duration), rational(tolerance)
    if duration <= 0 or tolerance < 0: raise ValueError('duration/tolerance')
    initial = clean(initial, gen.n, gen.max_terms)
    d = _denominator(gen.h, {'g': gen.gamma})
    q0 = _denominator(initial)
    if max(d.bit_length(), q0.bit_length()) > 8192: raise BudgetExceeded('denominator bit-length budget')
    b = {p: int(F(c) * q0) for p, c in initial.items() if c}
    coeffs = []
    probes = []
    gate_cost = dict(rejected_orders=0, inspected_max_entries=0, squared_entries=0, max_comparison_bits=0)
    sorting = pairs = conversions = 0
    columns = {}
    integer_cost = dict(multiply_adds=0, exported_fractions=0, max_integer_bits=0)
    def bounded(v):
        bits=abs(v).bit_length()
        integer_cost['max_integer_bits']=max(integer_cost['max_integer_bits'],bits)
        if bits>8192: raise BudgetExceeded('integer bit-length budget')
        return v
    for value in b.values(): bounded(value)
    # D action is d times the exact generator action, kept integral throughout.
    def apply_integer(op):
        nonlocal conversions
        out = {}
        for p, c in op.items():
            gen.cost['column_requests'] += 1
            if p not in columns:
                if len(columns) >= gen.max_terms: raise BudgetExceeded('cached column budget')
                gen.cost['hamiltonian_pairs'] += len(gen.h)
                col = commutator_i(gen.h, p)
                rate = -gen.gamma * sum(x != 'I' for x in p)
                if rate: col[p] = col.get(p, F(0)) + rate
                integer_col = {}
                for q, value in col.items():
                    z = value * d
                    if z.denominator != 1: raise ValueError('non-integral cleared column')
                    integer_col[q] = bounded(int(z))
                columns[p] = integer_col
            else: gen.cost['cache_hits'] += 1
            for q, z0 in columns[p].items():
                z = bounded(c * z0)
                out[q] = bounded(out.get(q, 0) + z)
                if not out[q]: del out[q]
                conversions += 1
                gen.cost['coefficient_multiply_adds'] += 1
                integer_cost['multiply_adds'] += 1
            if len(out) > gen.max_terms: raise BudgetExceeded('output term budget')
        gen.cost['peak_terms']=max(gen.cost['peak_terms'],len(op),len(out))
        return out
    for k in range(max_order + 1):
        denom = bounded(q0 * (d ** k) * factorial(k))
        c = {p: F(v, denom) for p, v in b.items() if v}
        coeffs.append(c)
        b_next = apply_integer(b)
        derivative_denom=bounded(denom*d)
        integer_cost['exported_fractions'] += len(c)
        gate = norm_rejection(b_next, derivative_denom, duration, tolerance, k, guard, gate_cost)
        if gate:
            probes.append(dict(order=k, grouping="rejection_guard", rejected_by=gate))
            gate_cost['rejected_orders'] += 1
            if k < max_order: b = b_next
            continue
        derivative = {p: F(v, derivative_denom) for p, v in b_next.items()}
        integer_cost['exported_fractions'] += len(derivative)
        for grouping in ('l1', 'firstfit', 'weighted'):
            groups, cost = norm_witness(derivative, grouping)
            sorting += cost['sorted_terms']; pairs += cost['group_comparisons']
            bound = duration ** (k + 1) / F(k + 1) * sum((F(g['upper']) for g in groups), F(0))
            probes.append(dict(order=k, grouping=grouping, bound=str(bound)))
            if bound <= tolerance:
                witness = {'schema':'v7-residual-1', 'integration_basis':'power',
                           'witnesses':{'jump:0':[], **{f'residual:0:{j}':[] for j in range(k + 1)},
                                        f'residual:0:{k}':groups if derivative else []},
                           'claimed_bound':str(bound),
                           'construction_cost':{'group_comparisons':pairs, 'sorted_terms':sorting,
                                                'integer_conversions':conversions, 'denominator':d, 'initial_denominator':q0}}
                if not derivative:
                    witness['witnesses']={'jump:0':[],'residual:0:0':[]}
                return Piece(duration, tuple(coeffs)), witness, dict(probes=probes, denominator=d, initial_denominator=q0, integer_work=integer_cost, generator=dict(gen.cost), guard_work=gate_cost)
        if k < max_order: b = b_next
    raise ValueError('Taylor order budget without certificate')


def norm_rejection(nums, denom, time, tolerance, order, guard, cost):
    """Necessary norm feasibility test. True means no valid upper bound can pass.

    Weight = time**(order+1)/(order+1), derivative = nums/denom.
    All comparisons use integers and strict >. Equality must proceed to the
    original norm proposals. Early exits and integer products are charged.
    """
    left = time.numerator**(order+1)*tolerance.denominator
    right = time.denominator**(order+1)*(order+1)*tolerance.numerator*denom
    bits=max(left.bit_length(),right.bit_length())
    cost['max_comparison_bits']=max(cost['max_comparison_bits'],bits)
    if bits>16384: raise BudgetExceeded('guard comparison bit budget')
    if guard in ('max','cascade'):
        threshold=right//left
        for value in nums.values():
            cost['inspected_max_entries']+=1
            if abs(value)>threshold:return 'max_coefficient'
    if guard in ('frobenius','cascade'):
        threshold=(right*right)//(left*left)
        total=0
        for value in nums.values():
            cost['squared_entries']+=1
            total+=value*value
            if total>threshold:return 'normalized_Hilbert_Schmidt'
    return None
