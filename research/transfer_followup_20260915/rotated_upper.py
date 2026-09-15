"""Exact rational orbital rotation + local-basis MPS upper, on the original H10.

The accepting rotation acts on one- and two-orbital operator coefficients,
never on the full fixed-N determinant space.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd, lcm
import os
from pathlib import Path
import subprocess
import time
from experiments.marginal_symbolic import decode, encode, canonical
from research.molecular_collective_20260913.core import digest
from research.transfer_solver_20260915.budget import ROOT, OUT, dump
from research.transfer_solver_20260915.campaign import STD, NUM, CHEM, PREFIX

SOURCE = OUT/'cases/h10_size'
CASE = OUT/'rotated_h10'


def mm(A, B):
    columns = list(zip(*B))
    return [[sum(a*b for a, b in zip(row, column) if a and b) for column in columns] for row in A]


def transpose(A):
    return [list(row) for row in zip(*A)]


def rational_orthogonal(target):
    import numpy as np
    n = len(target)
    A = target.copy()
    Z, den = [[int(i == j) for j in range(n)] for i in range(n)], 1
    for j in range(n-1):
        v = A[j:, j].copy()
        v[0] += np.copysign(np.linalg.norm(v), v[0])
        if max(abs(v)) < 1e-14:
            continue
        tail = [int(round(float(x/max(abs(v)))*10**8)) for x in v]
        vector = [0]*j+tail
        norm = sum(x*x for x in vector)
        P = [[(norm if a == b else 0)-2*vector[a]*vector[b] for b in range(n)] for a in range(n)]
        Z = mm(Z, P)
        den *= norm
        common = gcd(den, *(abs(x) for row in Z for x in row))
        den //= common
        Z = [[x//common for x in row] for row in Z]
        A = (np.array(P, dtype=float)/norm)@A
    signs = [1 if A[i, i] >= 0 else -1 for i in range(n)]
    Z = [[x*signs[j] for j, x in enumerate(row)] for row in Z]
    check_rotation(Z, den)
    return Z, den


def check_rotation(Z, den):
    n = len(Z)
    if type(den) is not int or den <= 0 or not n or any(len(row) != n or any(type(v) is not int for v in row) for row in Z):
        raise ValueError('Bounded exact square orbital rotation required')
    if den.bit_length() > 10000 or n > 12:
        raise ValueError('This orbital-rotation experiment is capped at 12 spatial orbitals')
    if mm(transpose(Z), Z) != [[den*den if i == j else 0 for j in range(n)] for i in range(n)]:
        raise ValueError('Orbital rotation is not exactly orthogonal')


def transformed(data, Z, den):
    check_rotation(Z, den)
    m = data['modes']
    if m != 2*len(Z):
        raise ValueError('Rotation does not match paired spin orbitals')
    h = decode(data['hamiltonian'], m, 4)
    hd = lcm(*(x.denominator for x in h.values()))
    rotation = [[Z[i//2][j//2] if i % 2 == j % 2 else 0 for j in range(m)] for i in range(m)]
    one = [[0]*m for _ in range(m)]
    pairs = list(combinations(range(m), 2))
    where = {pair: i for i, pair in enumerate(pairs)}
    two = [[0]*len(pairs) for _ in pairs]
    constant = 0
    for w, c in h.items():
        value = int(c*hd)
        if not w:
            constant += value
        elif len(w) == 2:
            if [kind for kind, mode in w] != [1, 0]:
                raise ValueError('Expected normal-ordered one-body input')
            one[w[0][1]][w[1][1]] += value
        else:
            if [kind for kind, mode in w] != [1, 1, 0, 0]:
                raise ValueError('Expected normal-ordered two-body input')
            p, q, s, r = [mode for kind, mode in w]
            sign = (1 if p < q else -1)*(1 if r < s else -1)
            two[where[tuple(sorted((p, q)))]][where[tuple(sorted((r, s)))]] += sign*value
    if one != transpose(one) or two != transpose(two):
        raise ValueError('One- and two-body coefficient matrices must be symmetric')
    transformed_one = mm(mm(transpose(rotation), one), rotation)
    # Paired spatial rotations preserve alpha count in each two-orbital pair.
    spin_groups = [[i for i, pair in enumerate(pairs) if sum(p % 2 == 0 for p in pair) == count] for count in range(3)]
    if any(two[i][j] for ga, group in enumerate(spin_groups) for gb, other in enumerate(spin_groups)
           if ga != gb for i in group for j in other):
        raise ValueError('This bounded transform requires alpha-count conservation')
    transformed_two = [[0]*len(pairs) for _ in pairs]
    for indices in spin_groups:
        W = [[rotation[p][i]*rotation[q][j]-rotation[p][j]*rotation[q][i]
              for i, j in (pairs[k] for k in indices)] for p, q in (pairs[k] for k in indices)]
        block = [[two[i][j] for j in indices] for i in indices]
        rotated = mm(mm(transpose(W), block), W)
        for a, i in enumerate(indices):
            for b, j in enumerate(indices):
                transformed_two[i][j] = rotated[a][b]
    return constant, transformed_one, transformed_two, hd, pairs


def rounded_division(numerator, denominator):
    quotient, remainder = divmod(numerator, denominator)
    return quotient+int(2*remainder > denominator or (2*remainder == denominator and quotient % 2))


def rotate_and_round(data, Z, den):
    constant, one, two, hd, pairs = transformed(data, Z, den)
    grid = 10**12
    common = hd*den**4
    error = 0
    terms = {}

    def accept(word, numerator):
        nonlocal error
        integer = rounded_division(numerator*grid, common)
        error += abs(numerator*grid-integer*common)
        if integer:
            terms[word] = terms.get(word, F(0))+F(integer, grid)

    accept((), constant*den**4)
    for i, row in enumerate(one):
        for j, value in enumerate(row):
            accept(((1, i), (0, j)), value*den**2)
    for a, row in enumerate(two):
        p, q = pairs[a]
        for b, value in enumerate(row):
            r, s = pairs[b]
            accept(((1, p), (1, q), (0, s), (0, r)), value)
    polynomial = canonical(terms)
    fixture = {'kind': 'rational_orthogonal_rotation_of_supplied_H_v1', 'modes': data['modes'],
        'particles': data['particles'], 'hamiltonian': encode(polynomial),
        'original_fixture_sha256': digest(data), 'orbital_rotation_denominator': str(den),
        'coefficient_rounding_norm_bound_Ha': str(F(error, common*grid)),
        'basis': 'Rational approximation to orthogonalized atomic orbitals; exact rotation of the supplied rational operator'}
    return fixture, F(error, common*grid)


def construct():
    import numpy as np
    from pyscf import gto
    started = time.monotonic()
    data = json.loads((SOURCE/'fixture.json').read_text())
    if (data['modes'], data['particles']) != (20, 10):
        raise ValueError('This size diagnostic is fixed to the supplied H10')
    z = np.load(SOURCE/'integrals.npz')
    mol = gto.M(atom=data['geometry'], basis=data['basis'], unit=data['unit'], verbose=0)
    overlap = mol.intor_symmetric('int1e_ovlp')
    eigenvalues, vectors = np.linalg.eigh(overlap)
    if min(eigenvalues) <= 1e-8:
        raise ValueError('Overlap is too ill conditioned for the declared orbital design')
    target = z['mo'].T@((vectors*np.sqrt(eigenvalues))@vectors.T)
    Z, den = rational_orthogonal(target)
    dump(CASE/'rotation.json', {'integer_matrix': Z, 'denominator': str(den),
        'original_fixture_sha256': digest(data), 'design': 'Symmetric AO orthogonalization in atom order',
        'floating_design_error': float(np.linalg.norm(np.array(Z, dtype=float)/den-target))})
    fixture, error = rotate_and_round(data, Z, den)
    dump(CASE/'fixture.json', fixture)
    dump(CASE/'rotation_construction.json', {'seconds': time.monotonic()-started, 'rounding_bound_Ha': str(error),
        'fixed_N_determinants_enumerated': 0, 'two_orbital_pair_coordinates': 190,
        'exact_orthogonality_checked': True, 'full_many_body_rotation_constructed': False})


def state():
    data = json.loads((CASE/'fixture.json').read_text())
    charges, tensors, current = [[[0, 0]]], [], [0, 0]
    for i in range(data['modes']):
        bit = int(i % 2 == (i//2) % 2)
        tensors.append([[0, bit, 0, 1]])
        current = current.copy()
        current[i % 2] += bit
        charges.append([current])
    seed = {'kind': 'integer_charge_mps_v1', 'fixture_sha256': digest(data), 'modes': 20, 'particles': 10,
        'spin_counts': [5, 5], 'denominator': 1, 'bond_charges': charges, 'tensors': tensors}
    dump(CASE/'neel_seed.json', seed)
    from research.correlated_pair_20260913.mps_spatial import run
    folder = CASE/'mps'
    folder.mkdir(exist_ok=False)
    run(data, [5, 5], folder, bond=32, sweeps=6, initial=str(CASE/'neel_seed.json'))


def verify():
    import sys
    from research.correlated_pair_20260913.mps_exact import check
    start = time.monotonic()
    data = json.loads((SOURCE/'fixture.json').read_text())
    rotation = json.loads((CASE/'rotation.json').read_text())
    if rotation['original_fixture_sha256'] != digest(data):
        raise ValueError('Original fixture binding failed')
    recomputed, error = rotate_and_round(data, rotation['integer_matrix'], int(rotation['denominator']))
    supplied = json.loads((CASE/'fixture.json').read_text())
    if recomputed != supplied:
        raise ValueError('The exact rotated operator and rounding allowance did not reproduce')
    path = CASE/'mps/state.json'
    if not path.exists():
        path = CASE/'mps/checkpoint.json'
    witness = json.loads(path.read_text())
    rec = check(supplied, witness)
    U = F(rec['upper_Ha'])+error
    if any(name in sys.modules for name in ('numpy', 'scipy', 'pyscf', 'quimb')):
        raise ValueError('Numerical import in orbital-rotation acceptance')
    dump(CASE/'original_upper.json', {'upper_Ha': str(U), 'upper_float_Ha': float(U),
        'original_fixture_sha256': digest(data), 'rotated_MPS_upper': rec,
        'rotation_rounding_allowance_Ha': str(error), 'exact_orbital_orthogonality': True,
        'used_state_path': str(path), 'full_fixed_N_space_enumerated': False,
        'two_particle_operator_coordinates': 190, 'numerical_imports': [], 'seconds': time.monotonic()-start,
        'argument': 'Second quantization of an exactly orthogonal paired orbital rotation is unitary and preserves N; the exact coefficient norm pays rounding.'})


def run():
    import hashlib
    CASE.mkdir(parents=True, exist_ok=False)
    steps = [('construct', 180, [CHEM, '-B', '-m', 'research.transfer_followup_20260915.rotated_upper', 'construct']),
             ('state', 600, [NUM, '-B', '-m', 'research.transfer_followup_20260915.rotated_upper', 'state']),
             ('verify', 600, [STD, '-B', '-S', '-m', 'research.transfer_followup_20260915.rotated_upper', 'verify'])]
    dump(CASE/'protocol.json', {'source_case': str(SOURCE), 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'old_many_body_state_used': False, 'FCI_used': False, 'steps': steps, 'bond_cap': 32, 'sweeps': 6,
        'purpose': 'Isolate a compact local-basis upper while keeping the original Hamiltonian and canonical lower family',
        'integral_generation_cost_is_shared_and_must_be_charged': True})
    env = os.environ.copy()
    env['NUMBA_CACHE_DIR'] = str(CASE/'numba_cache')
    outcomes = []
    for name, seconds, command in steps:
        result = subprocess.run([STD, '-B', '-S', '-m', PREFIX+'budget', '--name', 'h10_rotated_'+name,
                                '--seconds', str(seconds), '--']+command, cwd=ROOT, env=env)
        outcomes.append({'stage': name, 'exit_code': result.returncode})
        if result.returncode and name != 'state':
            break
        if result.returncode and name == 'state' and not (CASE/'mps/checkpoint.json').exists():
            break
    dump(CASE/'execution.json', outcomes)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('run', 'construct', 'state', 'verify'))
    a = p.parse_args()
    {'run': run, 'construct': construct, 'state': state, 'verify': verify}[a.action]()
