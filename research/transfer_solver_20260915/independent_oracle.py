"""Small, separately implemented exact bit-state validation of a full interval.

This is an explicitly enumerated validation control, not solver construction.
It imports no CAR/SOS or tensor-contraction checker implementation.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import time


def matrix(data):
    m, n = data['modes'], data['particles']
    states = [sum(1 << i for i in occupied) for occupied in combinations(range(m), n)]
    if len(states) > 100:
        raise ValueError('Independent exact validation oracle is limited to 100 determinants')
    position = {s: i for i, s in enumerate(states)}
    H = [[F(0) for _ in states] for _ in states]
    for j, s0 in enumerate(states):
        for term in data['hamiltonian']:
            s, sign = s0, 1
            for creation, orbital in reversed(term['word']):
                occupied = (s >> orbital) & 1
                if occupied == creation:
                    break
                if (s & ((1 << orbital)-1)).bit_count() % 2:
                    sign = -sign
                s ^= 1 << orbital
            else:
                H[position[s]][j] += sign*F(term['coefficient'])
    if any(H[i][j] != H[j][i] for i in range(len(states)) for j in range(i)):
        raise ValueError('Bit-state Hamiltonian is not exactly Hermitian')
    return states, H


def positive_shift(H, energy):
    n = len(H)
    L = [[F(0) for _ in range(n)] for _ in range(n)]
    D = []
    for i in range(n):
        pivot = H[i][i]-energy-sum(L[i][k]**2*D[k] for k in range(i))
        if pivot <= 0:
            return False, i
        D.append(pivot)
        L[i][i] = F(1)
        for j in range(i+1, n):
            L[j][i] = (H[j][i]-sum(L[j][k]*D[k]*L[i][k] for k in range(i)))/pivot
    return True, n


def amplitudes(states, state):
    answer = []
    for determinant in states:
        vector = {0: 1}
        for site, edges in enumerate(state['tensors']):
            bit = (determinant >> site) & 1
            following = {}
            for left, physical, right, coefficient in edges:
                if physical == bit and left in vector:
                    following[right] = following.get(right, 0)+vector[left]*coefficient
            vector = following
        answer.append(vector.get(0, 0))
    if not any(answer):
        raise ValueError('Zero independently expanded MPS')
    return answer


def run(case, output):
    start = time.monotonic()
    read = lambda p: json.loads(p.read_text())
    data, interval = read(case/'fixture.json'), read(case/'interval.json')
    states, H = matrix(data)
    L, U = F(interval['lower_Ha']), F(interval['upper_Ha'])
    accepted, pivots = positive_shift(H, L)
    if not accepted:
        raise ValueError('Independent exact LDL failed for claimed lower')
    vector = amplitudes(states, read(case/'mps/state.json'))
    rayleigh = sum(vector[i]*H[i][j]*vector[j] for i in range(len(states)) for j in range(len(states)))/sum(v*v for v in vector)
    if rayleigh != U:
        raise ValueError('Independent exact Rayleigh quotient differs from claimed upper')
    invalid_accepted, failed_pivot = positive_shift(H, U+F(1, 1000))
    if invalid_accepted:
        raise AssertionError('Invalid energy above a Rayleigh upper passed positivity')
    record = {'method': 'independent exact bit-action, rational LDL, and explicit integer MPS expansion',
              'full_sector_determinants_enumerated': len(states), 'positive_exact_pivots': pivots,
              'lower_Ha': str(L), 'upper_matches_exactly': True,
              'invalid_raised_lower_refused': True, 'failed_mutation_pivot': failed_pivot,
              'seconds': time.monotonic()-start,
              'enumeration_is_validation_only': True, 'original_CAR_or_MPS_checker_imported': False}
    with output.open('x') as stream:
        json.dump(record, stream, indent=2)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    run(a.case.resolve(), a.output.resolve())
