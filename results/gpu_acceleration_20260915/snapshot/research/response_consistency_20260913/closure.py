"""Exact low-dimensional obstruction to simple collective closure on H6."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import hashlib
import json
import sys
import time

from experiments.marginal_symbolic import decode
from research.molecular_collective_20260913.core import extract, factor_operators
from research.mechanism_transfer_20260913.ledger import modular_rank
from research.response_consistency_20260913.separator import act

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/response_consistency_20260913'


def primitive(matrix):
    den = lcm(*(F(x).denominator for row in matrix for x in row))
    integers = [[int(x*den) for x in row] for row in matrix]
    divisor = gcd(*(x for row in integers for x in row))
    return [[x//divisor for x in row] for row in integers] if divisor else integers


def commutator(a, b):
    n = len(a)
    return primitive([[sum(a[i][k]*b[k][j]-b[i][k]*a[k][j] for k in range(n))
                       for j in range(n)] for i in range(n)])


def lie_closure(generators):
    n = len(generators[0]); basis = []; provenance = []; rows = []
    if any(sum(a[i][i] for i in range(n)) for a in generators):
        raise ValueError('Traceless generators required for the stated upper dimension')
    attempts = 0
    def include(a, source):
        nonlocal attempts
        attempts += 1; row = [x for r in a for x in r]
        if modular_rank(rows+[row]) > len(rows):
            basis.append(a); rows.append(row); provenance.append(source)
    for i, a in enumerate(generators):
        include(primitive(a), {'generator': i})
    initial = len(basis); cursor = 0
    while cursor < len(basis) and len(basis) < n*n-1:
        for i, a in enumerate(generators):
            include(commutator(a, basis[cursor]), {'commutator': [i, cursor]})
            if len(basis) == n*n-1:
                break
        cursor += 1
    return {'initial_rank': initial, 'proved_rank_lower_bound': len(basis),
        'traceless_dimension_upper_bound': n*n-1, 'full_sl_closure_proved': len(basis) == n*n-1,
        'rank_trials': attempts, 'basis_provenance': provenance,
        'maximum_primitive_integer_bits': max(abs(x).bit_length() for row in rows for x in row),
        'method': 'Invertible rational reduction and pivots modulo 2^61-1 prove independence over Q. Matching the traceless upper dimension proves exact closure dimension.'}


def run():
    start = time.monotonic(); source = ROOT/'results/molecular_collective_20260913/campaign/h6'
    fixture_bytes = (source/'fixture.json').read_bytes(); tail_bytes = (source/'rank_10/tail.json').read_bytes()
    data = json.loads(fixture_bytes); tail = json.loads(tail_bytes); p = extract(data, tail['center_number'])
    patterns = [poly for _, poly in factor_operators(p, tail)]
    matrices = [[[poly.get(((1, 2*i), (0, 2*j)), F(0)) for j in range(6)] for i in range(6)] for poly in patterns]
    before = time.monotonic(); closure = lie_closure(matrices); closure['wall_seconds'] = time.monotonic()-before
    noncommuting = sum(any(x for row in commutator(a, b) for x in row) for a, b in combinations(matrices, 2))
    H = decode(data['hamiltonian'], 12, 4); occupations = ([0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 7]); diagonals = []
    for occupied in occupations:
        state = sum(1 << i for i in occupied)
        diagonals.append(act(H, {state: F(1)}).get(state, F(0)))
    gap = abs(diagonals[1]-diagonals[0])
    result = {'fixture_sha256': hashlib.sha256(fixture_bytes).hexdigest(), 'tail_sha256': hashlib.sha256(tail_bytes).hexdigest(),
        'spatial_pattern_count': len(patterns), 'noncommuting_pairs': noncommuting, 'pairs_checked': 45,
        'lie_closure': closure, 'occupation_split': {'eliminated_condition': 'spin orbital zero occupied',
            'diagnostic_occupations': occupations, 'exact_diagonals': [str(x) for x in diagonals],
            'diagonal_difference_Ha': str(gap), 'diagonal_difference_float_Ha': float(gap),
            'nonscalar_D_proved': bool(gap), 'many_body_determinants_evaluated': 2},
        'many_body_sector_matrices_enumerated': 0, 'wall_seconds': time.monotonic()-start,
        'scope': 'Refutes commuting-pattern closure and scalar D for this concrete occupancy split only. Does not exclude compact noncommutative response or prove a gap.'}
    if any(n in sys.modules for n in ('numpy', 'scipy', 'cvxpy', 'pyscf')):
        raise AssertionError('Numerical import during exact closure check')
    result['numerical_packages_loaded'] = []
    (OUT/'h6_closure.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
