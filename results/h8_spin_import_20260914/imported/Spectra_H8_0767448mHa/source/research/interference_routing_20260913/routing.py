"""Exact signed-network certificates; no numerical dependencies in replay.

The molecular adapter deliberately enumerates a bounded fixed-N sector. A
tree expresses all edge incidences as paths, and exact LDL checks their joint
signed capacity matrix. This is a finite oracle, not a compression theorem.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import argparse
import hashlib
import json
import time

from experiments.marginal_symbolic import decode, encode, hermitian
from research.all_angles_20260913.frontier.exact_h4_oracle import matrix, ldl_psd, rayleigh

MAX_STATES = 256
MAX_COMPONENT = 64


def rational(value):
    if type(value) not in (str, int, F):
        raise ValueError('Exact rational required')
    return F(value)


def load_sector(path):
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    m, n = data['modes'], data['particles']
    if type(m) is not int or type(n) is not int or not 0 <= n <= m <= 64:
        raise ValueError('Invalid physical sector')
    size = comb(m, n)
    if size > MAX_STATES:
        raise ValueError(f'Enumeration budget: {size} states exceeds {MAX_STATES}')
    h = decode(data['hamiltonian'], m, 4)
    if not hermitian(h) or any(sum(2*c-1 for c, _ in w) for w in h):
        raise ValueError('Hermitian number-conserving Hamiltonian required')
    a, states = matrix({'modes': m, 'hamiltonian': encode(h)}, n)
    return a, states, data, hashlib.sha256(raw).hexdigest()


def components(a):
    unseen = set(range(len(a)))
    parts = []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
        part = []
        while todo:
            i = todo.pop()
            part.append(i)
            for j in sorted(unseen):
                if a[i][j]:
                    unseen.remove(j)
                    todo.append(j)
        parts.append(sorted(part))
    return parts


def positive_tree(size, edges):
    """Deterministic maximum-weight positive spanning tree, for discovery."""
    parent = list(range(size))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    tree = []
    for i, j, w in sorted(edges, key=lambda e: (-e[2], e[0], e[1])):
        if w <= 0:
            continue
        x, y = root(i), root(j)
        if x != y:
            parent[x] = y
            tree.append([i, j])
    if len(tree) != size-1:
        raise ValueError('Positive subgraph disconnected')
    return tree


def path_capacity(size, edges, tree):
    """Return K in L = B_tree K B_tree^T, retaining every cross term."""
    if type(size) is not int or not 1 <= size <= MAX_COMPONENT+1:
        raise ValueError('Component work budget')
    weights = {}
    for i, j, raw in edges:
        if type(i) is not int or type(j) is not int or not 0 <= i < j < size:
            raise ValueError('Invalid edge endpoints')
        weights[i, j] = weights.get((i, j), F(0)) + rational(raw)
    if len(tree) != size-1:
        raise ValueError('Wrong tree size')
    adjacent = [[] for _ in range(size)]
    seen = set()
    for k, edge in enumerate(tree):
        if not isinstance(edge, list) or len(edge) != 2:
            raise ValueError('Invalid tree edge')
        i, j = edge
        if type(i) is not int or type(j) is not int or not 0 <= i < j < size:
            raise ValueError('Invalid tree endpoint')
        if (i, j) in seen or weights.get((i, j), F(0)) <= 0:
            raise ValueError('Tree must use distinct positive edges')
        seen.add((i, j))
        adjacent[i].append((j, k, 1))
        adjacent[j].append((i, k, -1))
    reachable = {0}
    todo = [0]
    while todo:
        for j, _, _ in adjacent[todo.pop()]:
            if j not in reachable:
                reachable.add(j)
                todo.append(j)
    if len(reachable) != size:
        raise ValueError('Disconnected tree')
    k_matrix = [[F(0)]*(size-1) for _ in range(size-1)]
    path_entries = 0
    products = 0
    for (start, end), w in weights.items():
        if not w:
            continue
        todo = [(start, -1, [])]
        route = None
        while todo:
            i, previous, path = todo.pop()
            if i == end:
                route = path
                break
            for j, k, sign in adjacent[i]:
                if j != previous:
                    todo.append((j, i, path+[(k, sign)]))
        if route is None:
            raise ValueError('Unroutable edge')
        path_entries += len(route)
        products += len(route)**2
        for i, x in route:
            for j, y in route:
                k_matrix[i][j] += w*x*y
    return k_matrix, {'edges': sum(bool(w) for w in weights.values()),
                      'negative_edges': sum(w < 0 for w in weights.values()),
                      'tree_edges': len(tree), 'path_entries': path_entries,
                      'capacity_additions': products}


def network(a, phi, lower=None):
    """With lower, add a ground vertex for the entire diagonal correction."""
    n = len(a)
    if len(phi) != n or any(type(x) is not int or x == 0 for x in phi):
        raise ValueError('Every lower-witness amplitude must be a nonzero integer')
    local = [sum(a[i][j]*phi[j] for j in range(n))/phi[i] for i in range(n)]
    edges = [(i, j, -a[i][j]*phi[i]*phi[j])
             for i in range(n) for j in range(i+1, n) if a[i][j]]
    if lower is not None:
        lower = rational(lower)
        edges += [(i, n, phi[i]**2*(local[i]-lower)) for i in range(n)
                  if local[i] != lower]
    return edges, local


def check_network(a, phi, tree, lower=None):
    edges, local = network(a, phi, lower)
    k_matrix, stats = path_capacity(len(a)+(lower is not None), edges, tree)
    ok, pivot, bad = ldl_psd(k_matrix)
    if not ok:
        raise ValueError(f'Joint signed capacity is not PSD; pivot {bad}')
    return min(local) if lower is None else rational(lower), {
        **stats, 'capacity_dimension': len(k_matrix),
        'minimum_ldl_pivot': str(pivot),
        'negative_local_slacks': 0 if lower is None else sum(x < rational(lower) for x in local)}


def replay(fixture, proof):
    start = time.monotonic()
    a, states, data, sha = load_sector(fixture)
    if proof['method'] not in ('split_tree_routing_v1', 'grounded_tree_routing_v1'):
        raise ValueError('Unknown proof method')
    if proof['fixture_sha256'] != sha or proof['modes'] != data['modes'] or proof['particles'] != data['particles']:
        raise ValueError('Fixture/sector mismatch')
    parts = components(a)
    if len(parts) != len(proof['components']):
        raise ValueError('Missing physical component')
    endpoints = []
    receipts = []
    for indices, witness in zip(parts, proof['components']):
        if len(indices) > MAX_COMPONENT:
            raise ValueError('Component work budget')
        if witness['states'] != [states[i] for i in indices]:
            raise ValueError('Wrong component states')
        block = [[a[i][j] for j in indices] for i in indices]
        lower = witness['lower'] if proof['method'] == 'grounded_tree_routing_v1' else None
        endpoint, stats = check_network(block, witness['amplitudes'], witness['tree'], lower)
        endpoints.append(endpoint)
        receipts.append(stats)
    v = proof['upper_amplitudes']
    if len(v) != len(a) or any(type(x) is not int for x in v) or not any(v):
        raise ValueError('Invalid variational upper witness')
    upper = rayleigh(a, v)
    lower = min(endpoints)
    if lower > upper:
        raise ValueError('Inconsistent interval')
    claims = {'lower': str(lower), 'upper': str(upper), 'width': str(upper-lower)}
    if 'claim' in proof and proof['claim'] != claims:
        raise ValueError('False stored interval')
    return {**claims, 'lower_float': float(lower), 'upper_float': float(upper),
            'width_float': float(upper-lower), 'passes_0_0016_Ha': upper-lower <= F(1, 625),
            'fixture_sha256': sha, 'dimension': len(a), 'component_sizes': list(map(len, parts)),
            'components': receipts, 'many_body_states_enumerated': len(a),
            'replay_seconds': time.monotonic()-start,
            'scope': 'Exact frozen finite Hamiltonian only; enumerative discovery and verification.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fixture', required=True, type=Path)
    parser.add_argument('--proof', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    result = replay(args.fixture, json.loads(args.proof.read_text()))
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'components'}))
