"""A four-mode exact obstruction to lifting a mode-graph PSD test."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import json

EDGES = ((0, 1, -1), (2, 3, -1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1))


def exchange_matrix(modes, particles, edges):
    states = [sum(1 << i for i in c) for c in combinations(range(modes), particles)]
    index = {s: i for i, s in enumerate(states)}
    a = [[F(0)]*len(states) for _ in states]
    for col, state in enumerate(states):
        for i, j, w in edges:
            if ((state >> i) ^ (state >> j)) & 1:
                a[col][col] += w
                a[index[state ^ (1 << i) ^ (1 << j)]][col] -= w
    return a, states


def run():
    one, _ = exchange_matrix(4, 1, EDGES)
    u = [1, 1, -1, -1]
    if one != [[F(x*y) for y in u] for x in u]:
        raise ValueError('One-particle outer-product identity failed')
    two, states = exchange_matrix(4, 2, EDGES)
    v = [{5: 1, 9: -1, 6: -1, 10: 1}.get(s, 0) for s in states]
    if any(sum(two[i][j]*v[j] for j in range(len(v))) != -2*v[i] for i in range(len(v))):
        raise ValueError('Negative two-particle eigenvector identity failed')
    return {'status': 'exact_counterexample', 'edges': EDGES,
            'one_particle_matrix': [[str(x) for x in row] for row in one],
            'one_particle_PSD_outer_factor': u, 'two_particle_states': states,
            'two_particle_matrix': [[str(x) for x in row] for row in two],
            'two_particle_integer_eigenvector': v, 'two_particle_eigenvalue': '-2',
            'scope': 'H=sum w_pq (I-SWAP_pq), equivalently parity-dressed CAR exchange. Not a counterexample to one-body second quantization dGamma(A).'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = run()
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'two_particle_eigenvalue': result['two_particle_eigenvalue']}))
