"""Bounded separator diagnostic using spin components of the same patterns.

This discovers a violated necessary positivity condition at the accepted dual
witness. It neither optimizes an energy nor proves a stronger cone optimal.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys
import time

from experiments.marginal_symbolic import add, mono, product, scale
from research.certificate_scaling.commutator_dual_witness import evaluate, moment_decode
from research.molecular_collective_20260913.core import digest, extract, factor_operators
from research.joint_patterns_20260913.core import anticommutator
from research.joint_patterns_20260913.dual import check as check_dual, cone_groups, integer_grams


def spin_generator(patterns, ref, modes):
    if not isinstance(ref, (list, tuple)) or len(ref) != 3 or any(type(v) is not int for v in ref):
        raise ValueError('Expected [pattern, spin, annihilation mode]')
    k, spin, mode = ref
    if not 0 <= mode < modes:
        raise ValueError('Invalid annihilation mode')
    a = mono(((0, mode),))
    if (k, spin) == (-1, -1):
        return a
    if not 0 <= k < len(patterns) or spin not in (0, 1):
        raise ValueError('Invalid pattern spin component')
    q = {w: c for w, c in patterns[k].items() if all(i % 2 == spin for _, i in w)}
    return product(q, a)


def check_separator(data, tail, dual, cert):
    start = time.monotonic()
    if cert.get('kind') != 'spin_pattern_separator_v1' or cert.get('fixture_sha256') != digest(data) or cert.get('tail_sha256') != digest(tail) or cert.get('dual_sha256') != digest(dual):
        raise ValueError('Spin separator input binding failed')
    p = extract(data, tail['center_number']); patterns = [q for _, q in factor_operators(p, tail)]
    refs = cert['generators']; vector = cert['vector']
    if not refs or len(refs) > p['modes']*(2*len(patterns)+1) or len(vector) != len(refs) or any(type(v) is not int for v in vector) or not any(vector):
        raise ValueError('Invalid bounded integer separator')
    q = add(*(scale(spin_generator(patterns, ref, p['modes']), c) for ref, c in zip(refs, vector) if c))
    poly = anticommutator(q, q)
    value = evaluate(poly, moment_decode(dual['moments'], p['modes']))
    if value >= 0:
        raise ValueError('This direction does not separate the dual witness')
    norm = sum(c*c for c in vector)
    return {'exact_expectation': str(value), 'coefficient_norm_squared': str(norm),
        'normalized_expectation': str(value/norm), 'normalized_expectation_float': float(value/norm),
        'generator_count': len(refs), 'nonzero_direction_coefficients': sum(c != 0 for c in vector),
        'operator_monomials': len(q), 'anticommutator_monomials': len(poly),
        'max_degree': max(map(len, poly), default=0), 'replay_seconds': time.monotonic()-start,
        'scope': 'Exact violation of a necessary positivity condition at this dual functional. No new energy bound or spin-resolved cone obstruction follows.'}


def propose(data, tail, dual, out):
    import numpy as np
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    accepted_dual = check_dual(data, tail, dual)
    p = extract(data, tail['center_number']); patterns = [q for _, q in factor_operators(p, tail)]
    _, signature = cone_groups(p, tail, dual['parity_masks']); parts = {}
    for k, spin in [(-1, -1)]+[(k, spin) for k in range(len(patterns)) for spin in range(2)]:
        for mode in range(p['modes']):
            ref = [k, spin, mode]; q = spin_generator(patterns, ref, p['modes'])
            if not q:
                continue
            keys = {signature(w) for w in q}
            if len(keys) != 1:
                raise ValueError('Spin-resolved generator lacks definite symmetry')
            group = parts.setdefault(next(iter(keys)), {'generators': [], 'polynomials': []})
            group['generators'].append(ref); group['polynomials'].append(q)
    y = moment_decode(dual['moments'], p['modes']); candidates = []; metrics = []
    for key in sorted(parts):
        group = parts[key]; matrices, den = integer_grams(group['polynomials'], [y], 'anticommutator')
        matrix = np.array([[float(F(v, den)) for v in row] for row in matrices[0]])
        values, vectors = np.linalg.eigh(matrix)
        metrics.append({'signature': key, 'dimension': len(values), 'minimum_proposed_eigenvalue': float(values[0])})
        candidates.append((values[0], group, np.rint(vectors[:, 0]*10**8).astype(np.int64)))
    _, group, direction = min(candidates, key=lambda item: item[0])
    cert = {'kind': 'spin_pattern_separator_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'dual_sha256': digest(dual), 'generators': group['generators'], 'vector': list(map(int, direction))}
    receipt = check_separator(data, tail, dual, cert)
    receipt.update({'proposal_blocks': metrics, 'wall_seconds': time.monotonic()-start,
        'dual_replay_seconds': accepted_dual['replay_seconds']})
    (out/'separator.json').write_text(json.dumps(cert, separators=(',', ':'))+'\n')
    (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt), flush=True); return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--fixture', type=Path, required=True)
    parser.add_argument('--tail', type=Path, required=True); parser.add_argument('--dual', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); parser.add_argument('--separator', type=Path)
    args = parser.parse_args(); data = json.loads(args.fixture.read_text()); tail = json.loads(args.tail.read_text()); dual = json.loads(args.dual.read_text())
    if args.separator:
        check_dual(data, tail, dual)
        result = check_separator(data, tail, dual, json.loads(args.separator.read_text()))
        if any(name in sys.modules for name in ('numpy', 'scipy', 'cvxpy', 'pyscf')):
            raise AssertionError('Numerical import on exact diagnostic replay')
        args.out.write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True)
    else:
        propose(data, tail, dual, args.out)
