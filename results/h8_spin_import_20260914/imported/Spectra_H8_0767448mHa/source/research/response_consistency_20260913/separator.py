"""Exact, existing-moment T1 separation; the old full audit stays separate."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json
import random
import sys
import time

from experiments.marginal_symbolic import add, canonical, mono, scale
from experiments.marginal_hunt_car import adj
from research.joint_patterns_20260913.core import anticommutator
from research.certificate_scaling.commutator_dual_witness import moment_decode, evaluate, seed

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/response_consistency_20260913'


def domain(m):
    return {tuple((1, i) for i in left)+tuple((0, i) for i in right)
            for k in range(3) for left in combinations(range(m), k)
            for right in combinations(range(m), k)}


def operator(cert):
    if cert.get('kind') != 'three_removal_separator_v1' or cert.get('modes') != 12:
        raise ValueError('This separator targets the twelve-mode parent')
    triples = cert['triples']; coefficients = cert['coefficients']; den = cert['denominator']
    if type(den) is not int or den <= 0 or not 1 <= len(triples) <= 220 or len(triples) != len(coefficients):
        raise ValueError('Invalid exact coefficient dimensions or denominator')
    seen = set()
    for triple, coefficient in zip(triples, coefficients):
        if len(triple) != 3 or any(type(i) is not int or not 0 <= i < 12 for i in triple):
            raise ValueError('Invalid annihilation triple')
        if sorted(set(triple)) != list(triple) or tuple(triple) in seen:
            raise ValueError('Triples must be distinct and increasing')
        if type(coefficient) is not int or not coefficient:
            raise ValueError('Nonzero integer coefficients required')
        seen.add(tuple(triple))
    if sorted({i//2 for triple in triples for i in triple}) != cert['support_spatial']:
        raise ValueError('Support does not describe the operator')
    C = add(*(scale(canonical(mono(tuple((0, i) for i in triple))), F(z, den))
              for triple, z in zip(triples, coefficients)))
    return C, anticommutator(C, C)


def check(cert, parent_bytes):
    start = time.monotonic()
    if hashlib.sha256(parent_bytes).hexdigest() != cert['parent_witness_sha256']:
        raise ValueError('Parent witness hash mismatch')
    parent = json.loads(parent_bytes)
    if parent.get('kind') != 'spin_completion_full_dual_v1':
        raise ValueError('Wrong parent family')
    C, P = operator(cert); words = domain(12); y = moment_decode(parent['moments'], 12)
    if any(w not in words for w in y) or any(w not in words for w in P):
        raise ValueError('Unspecified higher or unbalanced moments are forbidden')
    if y.get(()) != 1:
        raise ValueError('Parent is not normalized')
    value = evaluate(P, y)
    if value >= 0 or str(value) != cert['exact_moment']:
        raise ValueError('Claimed strict exact separation did not reproduce')
    return {'exact_moment': str(value), 'moment_float': float(value),
        'parent_witness_sha256': cert['parent_witness_sha256'],
        'operator_terms': len(C), 'expanded_polynomial_terms': len(P),
        'terms_by_degree': {str(k): sum(len(w) == k for w in P) for k in (0, 2, 4)},
        'specified_domain_words': len(words), 'sparse_zero_terms_used': sum(w not in y for w in P),
        'many_body_states_enumerated': 0, 'wall_seconds': time.monotonic()-start,
        'scope': 'Exact T1 separation of the hash-bound functional. Old-cone feasibility is checked independently by parent_family_audit.json.'}


def act(poly, vector):
    """Independent occupation-bit action for charged, small physical controls."""
    out = {}
    for state, amplitude in vector.items():
        for word, coefficient in poly.items():
            target = state; sign = 1
            for creation, mode in reversed(word):
                occupied = (target >> mode) & 1
                if occupied == creation:
                    break
                if (target & ((1 << mode)-1)).bit_count() % 2:
                    sign = -sign
                target ^= 1 << mode
            else:
                out[target] = out.get(target, F(0))+sign*coefficient*amplitude
    return {s: v for s, v in out.items() if v}


def controls(cert, parent):
    start = time.monotonic(); C, P = operator(cert); y = moment_decode(parent['moments'], 12)
    value = evaluate(P, y); trace_value = sum(c*seed(w, 12, 6) for w, c in P.items())
    threshold = -value/(trace_value-value); mixtures = []
    for t in (F(1, 1000), F(1, 200), F(1, 100), F(1, 50)):
        v = (1-t)*value+t*trace_value
        mixtures.append({'trace_fraction': str(t), 'exact_moment': str(v), 'rejected': v < 0})
    states = []
    for n in (4, 6, 8):
        basis = list(combinations(range(12), n)); rng = random.Random(20260913+n)
        vector = {sum(1 << i for i in inds): F(rng.choice((-3, -2, -1, 1, 2, 3)), 7)
                  for inds in rng.sample(basis, 32)}
        norm = sum(x*x for x in vector.values())
        left = act(C, vector); right = act(canonical(adj(C)), vector)
        squares = sum(x*x for x in left.values())+sum(x*x for x in right.values())
        expanded = sum(vector.get(s, F(0))*v for s, v in act(P, vector).items())
        if expanded != squares or squares < 0:
            raise AssertionError('Independent fermionic action disagrees with the positive identity')
        states.append({'particles': n, 'basis_states_enumerated': len(basis), 'trial_amplitudes': len(vector),
            'exact_expectation': str(expanded/norm), 'agrees_with_two_squared_norms': True})
    return {'uniform_fixed_N_trace_moment': str(trace_value),
        'strict_separation_for_trace_fraction_below': str(threshold),
        'threshold_float': float(threshold), 'fresh_controlled_mixtures': mixtures,
        'fresh_physical_controls': states, 'many_body_basis_states_enumerated': sum(x['basis_states_enumerated'] for x in states),
        'wall_seconds': time.monotonic()-start,
        'scope': 'Convex mixtures give an exact class of old-cone-feasible counterexamples. They are correlated controls, not independent molecular transfer.'}


def run():
    start = time.monotonic(); cert = json.loads((OUT/'separator.json').read_text())
    parent_bytes = (ROOT/'results/trace_pricing_20260913/full_dual/witness.json').read_bytes()
    result = check(cert, parent_bytes); result['controls'] = controls(cert, json.loads(parent_bytes))
    forbidden = [n for n in ('numpy', 'scipy', 'cvxpy', 'pyscf') if n in sys.modules]
    if forbidden:
        raise AssertionError('Numerical imports on the exact path')
    result['numerical_packages_loaded'] = forbidden; result['total_wall_seconds'] = time.monotonic()-start
    (OUT/'separator_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
