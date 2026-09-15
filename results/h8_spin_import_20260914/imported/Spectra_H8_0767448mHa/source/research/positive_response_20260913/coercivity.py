"""A rational, nonenumerating interaction-gap certificate from Lie brackets."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse
import json
import sys
import time

from research.molecular_collective_20260913.core import digest, extract, factor_operators, density, matmul
from research.mechanism_transfer_20260913.ledger import modular_rank
from experiments.marginal_symbolic import add, scale

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/positive_response_20260913'


def spatial_norm(a):
    return max(max(sum(abs(x) for x in row) for row in a),
               max(sum(abs(a[i][j]) for i in range(len(a))) for j in range(len(a))))


def coordinates(a):
    s = len(a)
    return [a[i][j] for i in range(s) for j in range(s) if i != j]+[a[i][i] for i in range(s-1)]


def spatial_inputs(data, tail):
    p = extract(data, tail['center_number']); s = p['spatial']; n = p['particles']
    if n != s:
        raise ValueError('Half filling required')
    matrices = []; weights = []
    for weight, q in factor_operators(p, tail):
        a = [[q.get(((1, 2*i), (0, 2*j)), F(0)) for j in range(s)] for i in range(s)]
        rebuilt = add(*(scale(density(i, j), a[i][j]) for i in range(s) for j in range(s) if a[i][j]))
        if q != rebuilt or weight <= 0 or sum(a[i][i] for i in range(s)):
            raise ValueError('Strictly positive, traceless spatial density squares required')
        matrices.append(a); weights.append(weight)
    return matrices, weights, n


def inputs(data, tail):
    result = spatial_inputs(data, tail)
    if result[2] != 6:
        raise ValueError('This quantitative certificate is capped at half-filled H6')
    return result


def node_value(node, basis, constants, matrices, weights, n):
    if set(node) == {'root'}:
        i = node['root']
        if type(i) is not int or not 0 <= i < len(matrices):
            raise ValueError('Invalid root pattern')
        raw = matrices[i]; r = spatial_norm(raw)
        if not r:
            raise ValueError('Zero root pattern')
        k = 2/(weights[i]*r*r)
    elif set(node) == {'bracket'}:
        pair = node['bracket']
        if len(pair) != 2 or any(type(i) is not int or not 0 <= i < len(basis) for i in pair):
            raise ValueError('Bracket must reference earlier nodes')
        i, j = pair; ab = matmul(basis[i], basis[j]); ba = matmul(basis[j], basis[i])
        raw = [[x-y for x, y in zip(row, other)] for row, other in zip(ab, ba)]
        r = spatial_norm(raw)
        if not r:
            raise ValueError('Zero bracket')
        k = 2*n*n*(constants[i]+constants[j])/(r*r)
    else:
        raise ValueError('Unknown Lie-basis node')
    a = [[x/r for x in row] for row in raw]
    if spatial_norm(a) != 1 or sum(a[i][i] for i in range(len(a))):
        raise AssertionError('Normalization or trace identity failed')
    return a, k


def inverse(a):
    size = len(a)
    work = [list(map(F, row))+[F(i == j) for j in range(size)] for i, row in enumerate(a)]
    for j in range(size):
        pivot = next((i for i in range(j, size) if work[i][j]), None)
        if pivot is None:
            raise ValueError('Singular Lie reconstruction')
        work[j], work[pivot] = work[pivot], work[j]; value = work[j][j]
        work[j] = [x/value for x in work[j]]
        for i in range(size):
            if i != j and work[i][j]:
                value = work[i][j]; work[i] = [x-value*y for x, y in zip(work[i], work[j])]
    return [row[size:] for row in work]


def sqrt_upper(x):
    den = 2**40; scaled = (x.numerator*den*den+x.denominator-1)//x.denominator
    root = isqrt(scaled)
    if root*root < scaled:
        root += 1
    return F(root, den)


def construct(data, tail):
    start = time.monotonic(); matrices, weights, n = inputs(data, tail)
    nodes = []; basis = []; constants = []; rows = []; root_ids = []
    attempts = 0
    def include(node):
        nonlocal attempts
        attempts += 1
        try:
            a, k = node_value(node, basis, constants, matrices, weights, n)
        except ValueError as error:
            if str(error) == 'Zero bracket':
                return False
            raise
        row = coordinates(a)
        if modular_rank(rows+[row]) > len(rows):
            nodes.append(node); basis.append(a); constants.append(k); rows.append(row); return True
        return False
    for i in range(len(matrices)):
        if include({'root': i}):
            root_ids.append(len(nodes)-1)
    cursor = 0
    while cursor < len(nodes) and len(nodes) < n*n-1:
        for i in root_ids:
            include({'bracket': [i, cursor]})
            if len(nodes) == n*n-1:
                break
        cursor += 1
    if len(nodes) != n*n-1:
        raise ValueError('Full traceless reconstruction was not established')
    frame = [list(row) for row in zip(*rows)]; inv = inverse(frame)
    cert = {'kind': 'lie_bulk_gap_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'nodes': nodes, 'sqrt_upper': [str(sqrt_upper(k)) for k in constants],
        'reconstruction': [[str(inv[j][i]) for j in range(len(inv))] for i in range(len(inv))]}
    stats = {'rank_trials': attempts, 'construction_seconds': time.monotonic()-start,
        'largest_matrix_dimension': len(nodes), 'many_body_states_enumerated': 0}
    return cert, stats


def check(data, tail, cert):
    start = time.monotonic(); matrices, weights, n = inputs(data, tail); size = n*n-1
    if cert.get('kind') != 'lie_bulk_gap_v1' or cert['fixture_sha256'] != digest(data) or cert['tail_sha256'] != digest(tail):
        raise ValueError('Gap certificate input binding failed')
    if len(cert['nodes']) != size or len(cert['sqrt_upper']) != size or len(cert['reconstruction']) != size:
        raise ValueError('Full bounded Lie-basis dimensions required')
    basis = []; constants = []
    for node in cert['nodes']:
        a, k = node_value(node, basis, constants, matrices, weights, n); basis.append(a); constants.append(k)
    if any(type(x) is not str for x in cert['sqrt_upper']):
        raise ValueError('Exact rational norm bounds required')
    upper = list(map(F, cert['sqrt_upper']))
    if any(u < 0 or u*u < k for u, k in zip(upper, constants)):
        raise ValueError('False propagated norm bound')
    coefficients = cert['reconstruction']
    if any(len(row) != size or any(type(x) is not str for x in row) for row in coefficients):
        raise ValueError('Exact reconstruction coefficients required')
    coefficients = [list(map(F, row)) for row in coefficients]; vectors = [coordinates(a) for a in basis]
    for target, row in enumerate(coefficients):
        if any(sum((c*v[j] for c, v in zip(row, vectors)), F(0)) != int(j == target) for j in range(size)):
            raise ValueError('Traceless target did not reconstruct exactly')
    K = sum(sum(abs(c)*u for c, u in zip(row, upper))**2 for row in coefficients)
    if K <= 0:
        raise ValueError('Positive coercivity constant required')
    return {'casimir_multiplier_lower': str(1/K), 'casimir_multiplier_lower_float': float(1/K),
        'double_occupancy_bulk_gap_Ha': str(2*n/K), 'double_occupancy_bulk_gap_float_Ha': float(2*n/K),
        'K': str(K), 'K_float': float(K), 'basis_dimension': size,
        'maximum_reconstruction_rational_bits': max(max(abs(c.numerator).bit_length(), c.denominator.bit_length()) for row in coefficients for c in row),
        'maximum_propagated_norm_squared_float': max(float(k) for k in constants),
        'replay_seconds': time.monotonic()-start, 'many_body_states_enumerated': 0,
        'scope': 'R>=C_orb/K on the fixed N=s sector, hence R>=2s/K on a specified doubly occupied orbital subspace. This is not yet a gap for the full molecular D.'}


def run(replay_only=False):
    OUT.mkdir(exist_ok=True); src = ROOT/'results/molecular_collective_20260913/campaign/h6'
    data = json.loads((src/'fixture.json').read_text()); tail = json.loads((src/'rank_10/tail.json').read_text())
    if replay_only:
        cert = json.loads((OUT/'gap_certificate.json').read_text()); receipt = check(data, tail, cert)
        receipt['numerical_packages_loaded'] = [x for x in ('numpy', 'scipy', 'cvxpy', 'pyscf') if x in sys.modules]
        if receipt['numerical_packages_loaded']:
            raise AssertionError('Numerical imports in exact gap replay')
        (OUT/'gap_replay.json').write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        cert, stats = construct(data, tail)
        text = json.dumps(cert, indent=2)+'\n'; (OUT/'gap_certificate.json').write_text(text)
        stats['certificate_bytes'] = len(text.encode()); (OUT/'gap_discovery.json').write_text(json.dumps(stats, indent=2)+'\n')
    print(json.dumps({'stage': 'replay' if replay_only else 'construction', 'complete': True}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
