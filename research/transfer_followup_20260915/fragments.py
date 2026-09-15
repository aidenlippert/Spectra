"""Exact composition for noninteracting fragments with fixed local charges.

This checks the tensor-sum rule. It is not a test of unrestricted charge
redistribution or of a monolithic SOS optimizer's size consistency.
"""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import decode, encode
from research.molecular_collective_20260913.core import digest
from research.correlated_pair_20260913.mps_exact import check as check_upper
from research.collective_completion_20260914.spin_screen import check as check_lower
from research.transfer_solver_20260915.budget import ROOT, OUT, dump


def construct(fragment, state, copies):
    m, n = fragment['modes'], fragment['particles']
    local = decode(fragment['hamiltonian'], m, 4)
    polynomial = {}
    for k in range(copies):
        for word, coefficient in local.items():
            shifted = tuple((creation, orbital+k*m) for creation, orbital in word)
            polynomial[shifted] = polynomial.get(shifted, F(0))+coefficient
    data = {'kind': 'noninteracting_fixed_fragment_charge_v1', 'modes': copies*m,
            'particles': copies*n, 'hamiltonian': encode(polynomial),
            'fragment_particle_constraints': [{'orbitals': list(range(k*m, (k+1)*m)), 'particles': n} for k in range(copies)],
            'interpretation': 'Fixed particle number on each fragment; not merely the fixed total-N space'}
    charges = [[[0, 0]]]
    for k in range(copies):
        for bond in state['bond_charges'][1:]:
            charges.append([[a+k*state['spin_counts'][0], b+k*state['spin_counts'][1]] for a, b in bond])
    combined = {**state, 'fixture_sha256': digest(data), 'modes': copies*m, 'particles': copies*n,
                'spin_counts': [copies*x for x in state['spin_counts']], 'bond_charges': charges,
                'tensors': state['tensors']*copies}
    return data, combined


def validate_decomposition(data, fragment, copies):
    m, n = fragment['modes'], fragment['particles']
    if (data['modes'], data['particles']) != (copies*m, copies*n):
        raise ValueError('Fragment sector dimensions do not add')
    expected = [{'orbitals': list(range(k*m, (k+1)*m)), 'particles': n} for k in range(copies)]
    if data.get('fragment_particle_constraints') != expected:
        raise ValueError('Local particle constraints are essential for this proof')
    local = decode(fragment['hamiltonian'], m, 4)
    wanted = {}
    for k in range(copies):
        for w, c in local.items():
            if sum(2*creation-1 for creation, orbital in w):
                raise ValueError('Each fragment must conserve particle number')
            shifted = tuple((creation, orbital+k*m) for creation, orbital in w)
            wanted[shifted] = wanted.get(shifted, F(0))+c
    if decode(data['hamiltonian'], data['modes'], 4) != wanted:
        raise ValueError('Hamiltonian does not equal the declared disjoint sum')


def run():
    started = time.monotonic()
    source = OUT/'cases/h4_control'
    read = lambda p: json.loads(p.read_text())
    fragment, state = read(source/'fixture.json'), read(source/'mps/state.json')
    upper = check_upper(fragment, state)
    lower = check_lower(fragment, read(source/'certificate.json'), read(source/'nonsinglet.json'))
    primitive_seconds = time.monotonic()-started
    U, L = F(upper['upper_Ha']), F(lower['lower'])
    directory = OUT/'fragments'
    directory.mkdir(exist_ok=False)
    results = []
    for copies in (2, 4, 8):
        t = time.monotonic()
        data, combined = construct(fragment, state, copies)
        validate_decomposition(data, fragment, copies)
        for k in range(1, copies+1):
            if combined['bond_charges'][k*fragment['modes']] != [[k*x for x in state['spin_counts']]]:
                raise ValueError('Upper state fails local charge boundary')
        replay = check_upper(data, combined)
        if F(replay['upper_Ha']) != copies*U:
            raise ValueError('Product-state exact upper failed additivity')
        cert = {'rule': 'Disjoint even number-conserving operator sum on fixed local-charge tensor factors',
                'primitive_fixture_sha256': digest(fragment), 'primitive_singlet_sha256': digest(read(source/'certificate.json')),
                'primitive_nonsinglet_sha256': digest(read(source/'nonsinglet.json')), 'copies': copies,
                'lower_Ha': str(copies*L), 'upper_Ha': str(copies*U),
                'width_mHa': float(1000*copies*(U-L)), 'product_upper_checked_directly': True,
                'local_lower_checked_once_and_reused': True, 'seconds_excluding_shared_primitive_replay': time.monotonic()-t}
        dump(directory/f'{copies}_fixture.json', data)
        dump(directory/f'{copies}_state.json', combined)
        dump(directory/f'{copies}_receipt.json', cert)
        results.append(cert)
    wrong = deepcopy(data)
    wrong.pop('fragment_particle_constraints')
    try:
        validate_decomposition(wrong, fragment, copies)
    except ValueError:
        pass
    else:
        raise AssertionError('Unrestricted total-charge reinterpretation was accepted')
    wrong = deepcopy(data)
    wrong['hamiltonian'].append({'word': [[1, 0], [0, fragment['modes']]], 'coefficient': '1/1000'})
    try:
        validate_decomposition(wrong, fragment, copies)
    except ValueError:
        pass
    else:
        raise AssertionError('Interfragment hopping was accepted')
    dump(directory/'summary.json', {'primitive_exact_replay_seconds': primitive_seconds,
        'complete_seconds': time.monotonic()-started, 'results': results,
        'wrong_local_sector_refused': True, 'interfragment_coupling_refused': True,
        'states_enumerated': 0, 'monolithic_family_size_consistency_tested': False,
        'claim': 'Certified intervals and product-state uppers compose additively for disjoint fragments with specified local particle numbers.'})


if __name__ == '__main__':
    run()
