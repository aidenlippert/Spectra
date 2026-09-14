"""Exact rational extension to finite-distance pair amplitudes on open chains.

This module verifies supplied witnesses. Numerical parameter discovery is separate.
The original nearest-neighbor implementation and its receipts remain unchanged.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time

from research.side_routes_20260913.positive_chain import rational, validate as validate_nearest


def validate(model, amplitude):
    m = model['modes']
    pairs = amplitude['pairs']
    if not isinstance(pairs, list) or len(pairs) > 3:
        raise ValueError('This experiment supports pair distances zero through three')
    if type(m) is not int or m < 2 or len(pairs) >= m:
        raise ValueError('Invalid chain size or pair range')
    base = {'sites': amplitude['sites'], 'bonds': pairs[0] if pairs else ['1'] * (m - 1)}
    m, n, t, v, h, y, _, offset = validate_nearest(model, base)
    xs = []
    for d, row in enumerate(pairs, 1):
        if not isinstance(row, list) or len(row) != m - d:
            raise ValueError('Pair-weight array length does not match its distance')
        values = [rational(a) for a in row]
        if any(a <= 0 for a in values):
            raise ValueError('Pair weights must be positive')
        xs.append(values)
    memory = min(m - 1, 2 * len(xs) + 1)
    if m * (n + 1) * (1 << memory) > 5000000:
        raise ValueError('Finite-range DP exceeds this runner state budget')
    return m, n, t, v, h, y, xs, offset, memory


def compile_factors(model, amplitude):
    m, n, t, v, h, y, xs, offset, memory = validate(model, amplitude)
    radius = len(xs)
    ending = [[(1, [F(0), field])] for field in h]
    entries = 2 * m
    for i in range(m - 1):
        start = max(0, i - radius)
        end = min(m - 1, i + 1 + radius)
        width = end - start + 1
        affected = []
        for d, weights in enumerate(xs, 1):
            for p in sorted({i - d, i, i + 1 - d, i + 1}):
                q = p + d
                if 0 <= p < q < m and {p, q} != {i, i + 1}:
                    affected.append((p, q, weights[p]))
        table = []
        for code in range(1 << width):
            bit = lambda j: (code >> (end - j)) & 1
            a, b = bit(i), bit(i + 1)
            value = v[i] * a * b
            if a != b:
                after = lambda j: b if j == i else a if j == i + 1 else bit(j)
                ratio = (y[i + 1] / y[i]) ** (a - b)
                for p, q, weight in affected:
                    ratio *= weight ** (after(p) * after(q) - bit(p) * bit(q))
                value -= t[i] * ratio
            table.append(value)
        ending[end].append((width, table))
        entries += len(table)
    return ending, entries


def bounds(model, amplitude):
    start = time.monotonic()
    m, n, t, v, h, y, xs, offset, memory = validate(model, amplitude)
    factors, entries = compile_factors(model, amplitude)
    states = {(0, 0): (offset, 0, F(1), offset)}
    tail_mask = (1 << memory) - 1
    transitions = 0
    peak_states = 1
    peak_bits = 1
    for i in range(m):
        nxt = {}
        remaining = m - i - 1
        for (count, tail), (minimum, witness, Z, S) in states.items():
            for bit in (0, 1):
                number = count + bit
                if number > n or number + remaining < n:
                    continue
                transitions += 1
                full = (tail << 1) | bit
                local = sum((table[full & ((1 << width) - 1)] for width, table in factors[i]), F(0))
                weight = y[i] ** (2 * bit)
                for d, weights in enumerate(xs, 1):
                    if i >= d:
                        weight *= weights[i - d] ** (2 * bit * ((tail >> (d - 1)) & 1))
                z, s = Z * weight, (S + Z * local) * weight
                lo, mask = minimum + local, (witness << 1) | bit
                key = (number, full & tail_mask)
                if key in nxt:
                    oldlo, oldmask, oldz, olds = nxt[key]
                    lo, mask = min((lo, mask), (oldlo, oldmask))
                    z += oldz
                    s += olds
                nxt[key] = lo, mask, z, s
                peak_bits = max(peak_bits, *(max(a.numerator.bit_length(), a.denominator.bit_length()) for a in (lo, z, s)))
        states = nxt
        peak_states = max(peak_states, len(states))
    lower, witness = min((a[0], a[1]) for a in states.values())
    Z = sum((a[2] for a in states.values()), F(0))
    S = sum((a[3] for a in states.values()), F(0))
    if Z <= 0:
        raise AssertionError('Nonpositive normalization')
    upper = S / Z
    if upper < lower:
        raise AssertionError('Reversed exact interval')
    return {'method': 'fixed_N_finite_pair_range_local_energy_DP_v1',
            'lower': str(lower), 'upper': str(upper), 'width': str(upper - lower),
            'lower_float': float(lower), 'upper_float': float(upper), 'width_float': float(upper - lower),
            'partition': str(Z), 'minimum_configuration': format(witness, f'0{m}b'),
            'modes': m, 'particles': n, 'pair_range': len(xs), 'tail_memory_bits': memory,
            'dp_transitions': transitions, 'peak_dp_states': peak_states,
            'factor_table_entries': entries, 'peak_dp_rational_bits': peak_bits,
            'many_body_states_enumerated': 0, 'wall_seconds': time.monotonic() - start,
            'scope': 'Exact supplied-chain interval; no global discovery optimum or generic accuracy guarantee.'}


def replay(payload):
    actual = bounds(payload['model'], payload['amplitude'])
    for key in ('lower', 'upper', 'width', 'partition', 'minimum_configuration'):
        if actual[key] != payload['claim'][key]:
            raise ValueError('False finite-range claim: ' + key)
    return actual


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--witness', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = replay(json.loads(args.witness.read_text()))
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
