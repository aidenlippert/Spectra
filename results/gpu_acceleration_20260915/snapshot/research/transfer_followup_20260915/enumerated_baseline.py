"""Small full-sector reference: exact rational residual-certified Cholesky.

All determinant labels and dense blocks are explicitly charged. This baseline
constructs its own upper and lower, with no Spectra state or lower as input.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb, lcm
from pathlib import Path
import time
from experiments.marginal_symbolic import decode
from research.certificate_scaling.streaming_reference_upper import compile_term
from research.compact_response_20260913.closure import check_factor
from research.molecular_collective_20260913.core import digest
from research.transfer_solver_20260915.budget import dump


def blocks(data):
    m, n = data['modes'], data['particles']
    if comb(m, n) > 2000:
        raise ValueError('Full determinant reference is bounded to 2000 labels')
    h = decode(data['hamiltonian'], m, 4)
    if any(sum(2*c-1 for c, i in w) or sum(2*c-1 for c, i in w if i % 2 == 0) for w in h):
        raise ValueError('Reference requires N and alpha count conservation')
    den = lcm(*(v.denominator for v in h.values()))
    terms = [term for w, c in h.items() if (term := compile_term(w, int(c*den))) is not None]
    groups = {}
    for occupied in combinations(range(m), n):
        state = sum(1 << i for i in occupied)
        groups.setdefault(sum(i % 2 == 0 for i in occupied), []).append(state)
    result = []
    for alpha, states in sorted(groups.items()):
        where = {s: i for i, s in enumerate(states)}
        H = [[0]*len(states) for _ in states]
        for j, state in enumerate(states):
            for required, occupied, flip, parity, value in terms:
                if state & required == occupied:
                    target = state ^ flip
                    if target not in where:
                        raise ValueError('Incomplete reference block')
                    H[where[target]][j] += value*(-1 if (state & parity).bit_count() % 2 else 1)
        if any(H[i][j] != H[j][i] for i in range(len(states)) for j in range(i)):
            raise ValueError('Reference block is not exactly Hermitian')
        result.append((alpha, states, H))
    return result, den, {'full_sector_labels_enumerated': comb(m, n),
                         'dense_matrix_entries': sum(len(s)**2 for a, s, H in result),
                         'largest_matrix_dimension': max(len(s) for a, s, H in result),
                         'Hamiltonian_word_state_checks': comb(m, n)*len(terms)}


def shifted(H, den, lower):
    kd = lcm(den, lower.denominator)
    return [[h*(kd//den)-(int(lower*kd) if i == j else 0) for j, h in enumerate(row)] for i, row in enumerate(H)], kd


def quotient(H, den, vector):
    if len(vector) != len(H) or any(type(v) is not int for v in vector) or not any(vector):
        raise ValueError('Nonzero complete integer reference vector required')
    return F(sum(vector[i]*h*vector[j] for i, row in enumerate(H) for j, h in enumerate(row)),
             den*sum(v*v for v in vector))


def construct(case):
    import numpy as np
    start = time.monotonic()
    data = json.loads((case/'fixture.json').read_text())
    groups, den, cost = blocks(data)
    rows = []
    for alpha, states, H in groups:
        hf = np.array(H, dtype=float)/den
        e, V = np.linalg.eigh(hf)
        lower = F(round(float(e[0])*10**10), 10**10)-F(8, 10000)
        fd = 10**12
        factor = np.linalg.cholesky(hf-(float(lower)+1e-7)*np.eye(len(H)))
        factor = [[int(round(float(factor[i, j])*fd)) for j in range(i+1)] for i in range(len(H))]
        vector = [int(round(float(x)*10**10)) for x in V[:, 0]]
        upper = quotient(H, den, vector)
        K, kd = shifted(H, den, lower)
        margin = check_factor(K, kd, F(0), factor, fd)
        rows.append({'alpha_count': alpha, 'lower_Ha': str(lower), 'upper_Ha': str(upper),
                     'factor_denominator': fd, 'factor': factor, 'upper_vector': vector,
                     'checked_margin_Ha': str(margin)})
    cert = {'kind': 'enumerated_reference_full_N_v1', 'fixture_sha256': digest(data), 'blocks': rows}
    directory = case/'enumerated_baseline'
    dump(directory/'certificate.json', cert)
    dump(directory/'construction.json', {'seconds': time.monotonic()-start, 'cost': cost,
        'existing_state_or_lower_used': False, 'target_width_Ha': '1/625',
        'planned_each_block_margin_Ha': '8/10000'})


def verify(case):
    import sys
    start = time.monotonic()
    data = json.loads((case/'fixture.json').read_text())
    directory = case/'enumerated_baseline'
    cert = json.loads((directory/'certificate.json').read_text())
    if cert.get('kind') != 'enumerated_reference_full_N_v1' or cert['fixture_sha256'] != digest(data):
        raise ValueError('Reference fixture binding failed')
    groups, den, cost = blocks(data)
    if [row['alpha_count'] for row in cert['blocks']] != [a for a, s, H in groups]:
        raise ValueError('Reference omitted a complete spin-projection block')
    lows, ups = [], []
    for (a, states, H), row in zip(groups, cert['blocks']):
        L = F(row['lower_Ha'])
        K, kd = shifted(H, den, L)
        margin = check_factor(K, kd, F(0), row['factor'], row['factor_denominator'])
        U = quotient(H, den, row['upper_vector'])
        if str(margin) != row['checked_margin_Ha'] or str(U) != row['upper_Ha'] or L > U:
            raise ValueError('Reference endpoints or remainder did not reproduce')
        lows.append(L)
        ups.append(U)
    L, U = min(lows), min(ups)
    if U-L > F(1, 625):
        raise ValueError('Reference missed the frozen total-width target')
    # A lower endpoint above a concrete Rayleigh quotient must be rejected.
    (a, states, H), row = next((g, r) for g, r in zip(groups, cert['blocks']) if F(r['upper_Ha']) == U)
    K, kd = shifted(H, den, U+F(1, 1000))
    try:
        check_factor(K, kd, F(0), row['factor'], row['factor_denominator'])
    except ValueError:
        pass
    else:
        raise AssertionError('Invalid raised lower was accepted')
    forbidden = [name for name in ('numpy', 'scipy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise ValueError('Numerical dependency in reference exact replay')
    dump(directory/'replay.json', {'lower_Ha': str(L), 'upper_Ha': str(U), 'width_mHa': float(1000*(U-L)),
        'target_met': True, 'cost': cost, 'seconds': time.monotonic()-start, 'valid_on': 'Entire fixed-N sector',
        'bad_lower_mutation_refused': True, 'numerical_imports': forbidden,
        'comparison_scope': 'Same rational fixture and <=1.6 mHa guarantee; determinant enumeration is explicit.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('construct', 'verify'))
    p.add_argument('case', type=Path)
    a = p.parse_args()
    construct(a.case.resolve()) if a.action == 'construct' else verify(a.case.resolve())
