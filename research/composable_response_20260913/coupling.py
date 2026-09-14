"""Local occupation blocks and exact CAR bounds; no determinant sector basis."""
from fractions import Fraction as F
from math import isqrt
import time

from experiments.marginal_symbolic import decode, canonical, add, word_product
from experiments.marginal_hunt_car import adj
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_general_schur import norm_bound


def bounded_product(p, q, deadline, stats):
    result = {}
    for left, a in p.items():
        for right, b in q.items():
            stats['word_products'] += 1
            if time.monotonic() > deadline:
                raise ValueError('Operator construction time cap exceeded')
            for word, sign in word_product(left, right):
                result[word] = result.get(word, F(0))+a*b*sign
            if len(result) > 100000:
                raise ValueError('Canonical operator word cap exceeded')
    return {w: c for w, c in result.items() if c}


def sqrt_up(value, denominator=10**8):
    if value < 0:
        raise ValueError('Negative squared norm')
    integer = isqrt(value.numerator*denominator**2//value.denominator)
    if F(integer**2, denominator**2) < value:
        integer += 1
    return F(integer, denominator)


def local_blocks(data):
    m = data['modes']
    if m not in (12, 16) or data['particles'] != m//2:
        raise ValueError('H6/H8 local occupation control required')
    r = m-4; inputs = [3, 7, 11]; outputs = [12, 13, 14, 15]
    polys = [[{} for _ in inputs] for _ in outputs]
    h = decode(data['hamiltonian'], m, 4)
    for word, coefficient in h.items():
        rest = tuple((c, i) for c, i in word if i < r)
        local = tuple((c, i-r) for c, i in word if i >= r)
        crossings = sum(i >= r and j < r for a, (_, i) in enumerate(word) for _, j in word[a+1:])
        for column, source in enumerate(inputs):
            value = apply_word(local, source)
            if value and value[0] in outputs:
                row = outputs.index(value[0]); p = polys[row][column]
                # Odd local words include the rest parity. Fixed total N fixes
                # this sign in each input local-occupation column.
                parity = crossings+len(local)*(data['particles']-source.bit_count())
                p[rest] = p.get(rest, F(0))+coefficient*value[1]*(-1 if parity % 2 else 1)
    return [[canonical(p) for p in row] for row in polys], inputs, outputs


def prove(data):
    start = time.monotonic(); deadline = start+60; stats = {'word_products': 0}
    blocks, inputs, outputs = local_blocks(data); norms = []; rows = []
    for row in blocks:
        nrow = []
        for p in row:
            if not p:
                nrow.append(F(0)); rows.append({'terms': 0, 'squared_norm_bound': '0'}); continue
            degrees = {len(w) % 2 for w in p}
            if len(degrees) != 1:
                raise ValueError('Mixed fermion parity in a local transition')
            square = bounded_product(adj(p), p, deadline, stats)
            if next(iter(degrees)):
                square = add(square, bounded_product(p, adj(p), deadline, stats))
            mu = norm_bound(square, 'hermitian_pairs'); nrow.append(sqrt_up(mu))
            rows.append({'terms': len(p), 'majorant_terms': len(square), 'squared_norm_bound': str(mu),
                'method': 'anticommutator' if next(iter(degrees)) else 'square'})
        norms.append(nrow)
    nu = max(map(sum, norms))*max(map(sum, zip(*norms)))
    return {'coupling_norm_squared_Ha2': str(nu), 'block_norm_bounds': [[str(x) for x in row] for row in norms],
        'local_input_occupations': inputs, 'local_output_occupations': outputs,
        'local_occupation_labels': 7, 'local_modes': 4, 'remaining_modes_not_enumerated': data['modes']-4,
        'blocks': rows, 'many_body_states_enumerated': 0, 'many_body_matrix_entries': 0,
        'word_products': stats['word_products'], 'wall_seconds': time.monotonic()-start}
