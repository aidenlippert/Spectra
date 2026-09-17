"""Exact product-operator to sparse MPO; no configurations or floating point.

A bipartite vertex cover factors each sparse coefficient matrix without
approximation or division. Input words are (creation_bit, mode).
"""
from collections import defaultdict, deque
from fractions import Fraction as F
import argparse
import json
import time
from pathlib import Path
from research.correlated_pair_20260913.mps_exact import local_word


def product_terms(fixture):
    m = fixture['modes']
    if type(m) is not int or m < 1:
        raise ValueError('Positive mode count required')
    out = defaultdict(F)
    for term in fixture['hamiltonian']:
        coeff = F(term['coefficient'])
        factors = []
        for mat in local_word(term['word'], m):
            first = next((x for x in mat if x), 0)
            if not first:
                coeff = F(0)
                break
            if first < 0:
                mat = tuple(-x for x in mat)
                coeff = -coeff
            factors.append(mat)
        if coeff:
            out[tuple(factors)] += coeff
    return {w: c for w, c in out.items() if c}


def vertex_cover(adjacency):
    match_right = {}
    def augment(left, seen):
        for right in sorted(adjacency[left]):
            if right in seen:
                continue
            seen.add(right)
            if right not in match_right or augment(match_right[right], seen):
                match_right[right] = left
                return True
        return False
    for left in sorted(adjacency):
        augment(left, set())
    match_left = {left: right for right, left in match_right.items()}
    reached_left = set(adjacency) - set(match_left)
    reached_right = set()
    queue = deque(sorted(reached_left))
    while queue:
        left = queue.popleft()
        for right in adjacency[left]:
            if match_left.get(left) == right or right in reached_right:
                continue
            reached_right.add(right)
            if right in match_right and match_right[right] not in reached_left:
                reached_left.add(match_right[right])
                queue.append(match_right[right])
    return set(adjacency) - reached_left, reached_right


def from_products(terms, modes):
    if any(len(w) != modes for w in terms):
        raise ValueError('Product operator length')
    states = {(0, w): F(c) for w, c in terms.items() if c}
    layers = []
    widths = [1]
    for site in range(modes):
        edges = defaultdict(F)
        for (left, word), coeff in states.items():
            edges[((left, word[0]), word[1:])] += coeff
        edges = {key: c for key, c in edges.items() if c}
        layer = defaultdict(F)
        following = defaultdict(F)
        if site == modes - 1:
            for ((left, mat), suffix), coeff in edges.items():
                assert suffix == ()
                layer[left, 0, mat] += coeff
            width = 1
        else:
            adjacency = defaultdict(set)
            for row, suffix in edges:
                adjacency[row].add(suffix)
            left_cover, right_cover = vertex_cover(adjacency)
            row_ids = {row: i for i, row in enumerate(sorted(left_cover))}
            for (left, mat), idx in row_ids.items():
                layer[left, idx, mat] = F(1)
            used_columns = sorted({suffix for row, suffix in edges if row not in left_cover})
            assert set(used_columns) <= right_cover
            col_ids = {suffix: len(row_ids) + i for i, suffix in enumerate(used_columns)}
            for suffix, idx in col_ids.items():
                following[idx, suffix] = F(1)
            for (row, suffix), coeff in edges.items():
                if row in row_ids:
                    following[row_ids[row], suffix] += coeff
                else:
                    left, mat = row
                    layer[left, col_ids[suffix], mat] += coeff
            width = len(row_ids) + len(col_ids)
        layers.append([[left, right, list(mat), str(c)]
                       for (left, right, mat), c in sorted(layer.items()) if c])
        widths.append(width)
        states = {key: c for key, c in following.items() if c}
    return {'kind': 'exact_real_product_mpo_v1', 'modes': modes,
            'widths': widths, 'layers': layers,
            'product_terms': len(terms), 'enumerated_configurations': 0}


def build(fixture):
    return from_products(product_terms(fixture), fixture['modes'])


def expand_products(mpo, max_terms=100000):
    """Symbolic operator reconstruction; no physical state amplitudes."""
    states = {(0, ()): F(1)}
    for layer in mpo['layers']:
        rows = defaultdict(list)
        for left, right, mat, c in layer:
            rows[left].append((right, tuple(mat), F(c)))
        following = defaultdict(F)
        for (left, word), value in states.items():
            for right, mat, c in rows[left]:
                following[right, word + (mat,)] += value * c
        states = {key: c for key, c in following.items() if c}
        if len(states) > max_terms:
            raise ValueError('Symbolic reconstruction term budget exceeded')
    if any(right != 0 for right, _ in states):
        raise ValueError('Non-scalar MPO boundary')
    return {word: c for (_, word), c in states.items()}


def sparse_mpo_edges(mpo):
    """Layers of (left,right,bra,ket,Fraction coefficient)."""
    layers = []
    for layer in mpo['layers']:
        merged = defaultdict(F)
        for left, right, mat, c in layer:
            for k, value in enumerate(mat):
                if value:
                    merged[left, right, k // 2, k % 2] += F(c) * value
        layers.append([(*key, value) for key, value in sorted(merged.items()) if value])
    return layers


def element(mpo, bra, ket):
    env = {0: F(1)}
    for site, layer in enumerate(mpo['layers']):
        out = defaultdict(F)
        bit = 2 * ((bra >> site) & 1) + ((ket >> site) & 1)
        for left, right, mat, c in layer:
            if left in env and mat[bit]:
                out[right] += env[left] * F(c) * mat[bit]
        env = {key: value for key, value in out.items() if value}
    return env.get(0, F(0))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('fixture'); p.add_argument('--output', required=True)
    args = p.parse_args()
    start = time.monotonic()
    data = json.loads(Path(args.fixture).read_text())
    op = build(data)
    if expand_products(op) != product_terms(data):
        raise AssertionError('Exact symbolic reconstruction failed')
    op['construction_and_symbolic_replay_seconds'] = time.monotonic() - start
    out = Path(args.output)
    if out.exists():
        raise ValueError('Refuse to overwrite an existing result')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(op, separators=(',', ':')) + '\n')
    print(json.dumps({key: value for key, value in op.items() if key != 'layers'}))
