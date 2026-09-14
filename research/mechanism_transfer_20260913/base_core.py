"""Exact learned-pattern anticommutators and transfer to a frozen molecule.

Accepting paths use only the standard library and the existing rational CAR
verifier. A compact factor is expanded once and its adjoint is constructed
exactly: rounding cannot untie the two halves of an anticommutator.
"""
from fractions import Fraction as F
from math import lcm
import json
import time

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import (
    add, canonical, decode, encode, mono, product, scale, validate_word, verify,
    word_product,
)
from research.molecular_collective_20260913.core import (
    digest, extract, factor_operators, retained_polynomial, tail_replay,
)


def anticommutator(left, right):
    """{left^dagger,right}; cubic odd leading terms must cancel, never drop."""
    for p in (left, right):
        if not p or canonical(p) != p or any(len(w) not in (1, 3) for w in p):
            raise ValueError('Canonical odd polynomials of degree at most three required')
    dagger = canonical(adj(left))
    out = add(product(dagger, right), product(right, dagger))
    if any(len(w) > 4 for w in out):
        raise AssertionError('Six-operator anticommutator cancellation failed')
    return out


def prepare_anticommutators(groups):
    """Exact monomial map in the existing contraction backend's format."""
    blocks = []; allwords = set(); pairs = 0; nonzeros = 0
    for group in groups:
        words = sorted({w for p in group['polynomials'] for w in p}, key=lambda w: (len(w), w))
        if not words or any(len(w) not in (1, 3) for w in words):
            raise ValueError('Odd generator words required')
        size = len(words); pairs += size * size; terms = []
        for i, left in enumerate(words):
            dagger = tuple((1-c, j) for c, j in reversed(left))
            for j, right in enumerate(words):
                reduced = {}
                for a, b in ((dagger, right), (right, dagger)):
                    for w, c in word_product(a, b):
                        reduced[w] = reduced.get(w, 0) + c
                for w, c in reduced.items():
                    if not c:
                        continue
                    if len(w) > 4:
                        raise AssertionError('Six-operator terms did not cancel exactly')
                    terms.append((w, i*size+j, c)); allwords.add(w)
        blocks.append({'words': words, 'terms': terms}); nonzeros += len(terms)
    return blocks, allwords, {'monomial_word_pairs': pairs,
        'CAR_word_products': 2*pairs, 'monomial_map_nonzeros': nonzeros,
        'max_surviving_degree': max(map(len, allwords), default=0)}


def generator(patterns, ref, modes):
    if not isinstance(ref, (list, tuple)) or len(ref) != 2:
        raise ValueError('A generator reference is [pattern, annihilation mode]')
    k, p = ref
    if type(k) is not int or type(p) is not int or not -1 <= k < len(patterns) or not 0 <= p < modes:
        raise ValueError('Invalid learned-pattern generator reference')
    a = mono(((0, p),))
    return a if k == -1 else product(patterns[k], a)


def pack_rows(name, rows, denominator):
    words = sorted({w for p in rows for w in p}, key=lambda w: (len(w), w))
    factors = []
    for p in rows:
        row = [p.get(w, F(0))*denominator for w in words]
        if any(F(c).denominator != 1 for c in row):
            raise ValueError('Common denominator does not represent the factors exactly')
        factors.append([int(c) for c in row])
    return {'name': name, 'words': words, 'factor': factors}


def expand_certificate(data, tail, cert):
    if cert.get('kind') != 'joint_density_anticommutator_v1':
        raise ValueError('Unknown joint-pattern certificate')
    if cert.get('fixture_sha256') != digest(data) or cert.get('tail_sha256') != digest(tail):
        raise ValueError('Frozen molecule or retained-pattern binding failed')
    m = data['modes']
    if type(m) is not int or not 4 <= m <= 16:
        raise ValueError('This transfer pass is capped at sixteen spin orbitals')
    den = cert['denominator']
    if type(den) is not int or den <= 0:
        raise ValueError('Positive integer factor denominator required')
    decode(cert['number_multiplier'], m, 2)
    p = extract(data, tail['center_number'])
    patterns = [a for _, a in factor_operators(p, tail)]
    expanded = []; common = den
    for block in cert['base_blocks']:
        # The baseline is quadratic. The ordinary verifier handles charge,
        # integer factors, and all remaining input validation after expansion.
        words = [validate_word(w, m, 2) for w in block['words']]
        expanded.append((block['name'], None, {**block, 'words': words}))
    if len(cert['anti_blocks']) > 4*(len(patterns)+1):
        raise ValueError('Anticommutator block count exceeds this pass budget')
    for block in cert['anti_blocks']:
        refs = block['generators']
        if not refs or len(refs) > m*(len(patterns)+1):
            raise ValueError('Invalid anticommutator generator count')
        polys = [generator(patterns, ref, m) for ref in refs]
        factors = block['factor']
        if len(factors) > len(refs) or any(len(row) != len(refs) or any(type(c) is not int for c in row) for row in factors):
            raise ValueError('Invalid compact integer factor')
        rows = []
        for row in factors:
            q = add(*(scale(poly, F(c, den)) for c, poly in zip(row, polys) if c))
            if q:
                rows.append(q)
                for c in q.values():
                    common = lcm(common, F(c).denominator)
        if rows:
            expanded.append((block['name']+'/B', rows, None))
            expanded.append((block['name']+'/B_dagger', [canonical(adj(q)) for q in rows], None))
    blocks = []
    for name, rows, block in expanded:
        if rows is None:
            if any(len(row) != len(block['words']) or any(type(c) is not int for c in row) for row in block['factor']):
                raise ValueError('Invalid baseline integer factor')
            blocks.append({**block, 'factor': [[c*(common//den) for c in row] for row in block['factor']]})
        else:
            blocks.append(pack_rows(name, rows, common))
    return {'modes': m, 'particles': data['particles'],
        'hamiltonian': encode(retained_polynomial(p, tail)),
        'b': cert['b'], 'number_multiplier': cert['number_multiplier'],
        'denominator': common, 'blocks': blocks}


def replay(data, tail, cert):
    start = time.monotonic()
    tail_receipt = tail_replay(data, tail)
    expanded = expand_certificate(data, tail, cert)
    built = time.monotonic()
    exact = verify(expanded)
    if exact['residual_max_degree'] > 4:
        raise AssertionError('Paired exact factors left a degree-six residual')
    original_lower = F(exact['lower']) + F(tail_receipt['lower_operator_shift_Ha'])
    return {'retained': exact, 'tail': tail_receipt,
        'original_lower_Ha': str(original_lower), 'original_lower_float_Ha': float(original_lower),
        'expand_seconds': built-start, 'SOS_replay_seconds': time.monotonic()-built,
        'replay_seconds': time.monotonic()-start,
        'expanded_certificate_bytes': len(json.dumps(expanded, separators=(',', ':')).encode())+1,
        'many_body_states_enumerated': 0,
        'scope': 'Exact lower bound on the frozen molecular Hamiltonian in its fixed-N sector. No cone-optimality claim.'}
