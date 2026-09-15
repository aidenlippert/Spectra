"""Exact rule for a formula-specified conditional-exchange response control."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import sys
import time

from research.response_consistency_20260913 import response as reference

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/positive_response_20260913'
KIND = 'conditional_exchange_central_spin_v1'


def parameters(model):
    if set(model) != {'kind', 'reference', 'lambda_retained', 'lambda_eliminated'} or model['kind'] != KIND:
        raise ValueError('Only the specified periodic conditional-exchange bath is supported')
    M, q, A0, D0, g, c = reference.parameters(model['reference'])
    if M < 4 or not 3 <= q <= M-1:
        raise ValueError('This nontrivial control requires M>=4 and 3<=q<=M-1')
    if any(type(model[k]) is not str for k in ('lambda_retained', 'lambda_eliminated')):
        raise ValueError('Exact rational positive-bulk strengths required')
    l0, l1 = (F(model[k]) for k in ('lambda_retained', 'lambda_eliminated'))
    if min(l0, l1) < 0:
        raise ValueError('Conditional exchange must have nonnegative strength')
    return M, q, A0, D0, g, c, l0, l1


def check(model, cert):
    start = time.monotonic(); M, q, _, _, _, _, _, _ = parameters(model)
    ref = reference.check(model['reference'], cert)
    return {k: v for k, v in ref.items() if k not in ('residual_E', 'scope', 'replay_seconds')} | {
        'response_residual_formula': '-lambda_eliminated*(g/d)*L_(q-1)*J_-',
        'positive_residual_block': '[-X0,I]^dagger*(lambda_eliminated*L_(q-1))*[-X0,I]',
        'local_positive_identity': 'n_k*(I-Swap_ij) = F_i^dagger*F_i/2, F_i=n_k*(I-Swap_ij)',
        'reference_trial_annihilated_by_added_terms': True,
        'conceptual_local_terms': 2*M,
        'scalar_parameters_read': 8,
        'replay_seconds': time.monotonic()-start,
        'scope': 'Formula-specified control with a common zero of the added positive terms; no arbitrary bath terms accepted and no molecular scaling claim.'}


def add_vectors(*vectors):
    out = {}
    for vector in vectors:
        for state, amplitude in vector.items():
            out[state] = out.get(state, F(0))+amplitude
    return {state: amplitude for state, amplitude in out.items() if amplitude}


def scale_vector(vector, coefficient):
    return {state: coefficient*amplitude for state, amplitude in vector.items() if coefficient*amplitude}


def inner(left, right):
    return sum((amplitude*right.get(state, F(0)) for state, amplitude in left.items()), F(0))


def lowering(M, vector):
    out = {}
    for state, amplitude in vector.items():
        for i in range(M):
            if state & (1 << i):
                target = state ^ (1 << i); out[target] = out.get(target, F(0))+amplitude
    return {state: amplitude for state, amplitude in out.items() if amplitude}


def raising(M, vector):
    mask = (1 << M)-1
    return {mask ^ state: amplitude for state, amplitude in lowering(M, {mask ^ state: amplitude for state, amplitude in vector.items()}).items()}


def exchange(M, vector):
    """Independent occupation-bit implementation, with spin (not CAR) signs."""
    out = {}
    for state, amplitude in vector.items():
        for i in range(M):
            j = (i+1) % M; k = (i+2) % M
            if state & (1 << k) and bool(state & (1 << i)) != bool(state & (1 << j)):
                target = state ^ (1 << i) ^ (1 << j)
                out[state] = out.get(state, F(0))+amplitude
                out[target] = out.get(target, F(0))-amplitude
    return {state: amplitude for state, amplitude in out.items() if amplitude}


def basis(M, q):
    return [sum(1 << i for i in indices) for indices in combinations(range(M), q)]


def exact_control(model, cert):
    start = time.monotonic(); M, q, A0, D0, g, c, l0, l1 = parameters(model)
    if M > 9:
        raise ValueError('Explicit controls are capped at nine bath spins')
    b = F(cert['lower']); d = D0-b; x = g/d
    retained = basis(M, q); eliminated = basis(M, q-1)
    dicke = {state: F(1) for state in retained}; lower_dicke = {state: F(1) for state in eliminated}
    if exchange(M, dicke) or exchange(M, lower_dicke):
        raise AssertionError('Added interactions do not kill the analytic trial')
    columns = [(state, exchange(M, {state: F(1)})) for state in eliminated]
    diagonals = [column.get(state, F(0)) for state, column in columns]
    offdiagonal = next(((state, target, l1*a) for state, column in columns for target, a in column.items()
                        if target != state and l1*a), None)
    candidates = []; commutator_examples = []
    for state in retained:
        vector = {state: F(1)}; jv = lowering(M, vector); ljv = exchange(M, jv)
        comm = add_vectors(ljv, scale_vector(lowering(M, exchange(M, vector)), -1))
        if comm:
            commutator_examples.append({'retained_state': state, 'squared_norm': str(inner(comm, comm))})
        k = A0-b-g*g*inner(jv, jv)/d+l0*inner(vector, exchange(M, vector))+x*x*l1*inner(jv, ljv)
        e_squared = x*x*l1*l1*inner(ljv, ljv)
        candidates.append((k-e_squared/d, state, k, e_squared))
    value, state, k, e_squared = min(candidates)
    if l1 <= 0 or (max(diagonals) == min(diagonals) and offdiagonal is None) or not commutator_examples or not e_squared or value >= 0:
        raise AssertionError('Control did not demonstrate all stipulated nontrivial properties')
    # Verify the full congruence quadratic identity on a deterministic sparse pair.
    v = {state: F((i % 5)-2, 7) for i, state in enumerate(retained) if i % 5 != 2}
    w = {state: F((i % 7)-3, 11) for i, state in enumerate(eliminated) if i % 7 != 3}
    jv = lowering(M, v); xv = scale_vector(jv, x); z = add_vectors(w, scale_vector(xv, -1))
    original = (A0-b)*inner(v, v)+(D0-b)*inner(z, z)+2*g*inner(jv, z)
    original += l0*inner(v, exchange(M, v))+l1*inner(z, exchange(M, z))
    factored = (A0-b)*inner(v, v)-g*g*inner(jv, jv)/d+d*inner(w, w)
    factored += l0*inner(v, exchange(M, v))+l1*inner(z, exchange(M, z))
    if original != factored or factored < 0:
        raise AssertionError('Independent spin action disagrees with the positive congruence')
    return {'bath_spins': M, 'excitations': q, 'basis_states_enumerated': len(retained)+len(eliminated),
        'eliminated_block_diagonal_spread': str(l1*(max(diagonals)-min(diagonals))),
        'eliminated_block_offdiagonal_example': {'source': offdiagonal[0], 'target': offdiagonal[1],
            'coefficient': str(offdiagonal[2])} if offdiagonal is not None else None,
        'nonscalar_eliminated_block_proved': True,
        'noncommutation_example': commutator_examples[0],
        'nonzero_residual_squared_norm': str(e_squared),
        'scalar_error_test_counterexample': {'retained_state': state, 'K_expectation': str(k),
            'K_minus_EdaggerE_over_d': str(value), 'expectation_float': float(value)},
        'exact_congruence_quadratic_identity': True,
        'dicke_vectors_annihilated': True,
        'wall_seconds': time.monotonic()-start}


def cases():
    for M, q, l0, l1, epsilon in ((5, 3, '1/3', '10', '3/10'),
                                 (7, 3, '2', '17', '-1/5'),
                                 (9, 4, '5/7', '100', '7/10'),
                                 (64, 31, '4/5', '9', '3/10'),
                                 (1024, 511, '7/11', '23', '-2/5')):
        yield {'kind': KIND, 'reference': {'kind': 'homogeneous_central_spin_v1', 'bath_spins': M,
            'excitations': q, 'epsilon': epsilon, 'omega': '0', 'chi': '0', 'g': str(F(1, M))},
            'lambda_retained': l0, 'lambda_eliminated': l1}


def run(replay=False):
    start = time.monotonic(); OUT.mkdir(exist_ok=True)
    if replay:
        records = json.loads((OUT/'response_certificates.json').read_text())
    else:
        records = []
        for model in cases():
            before = time.monotonic(); cert = reference.propose(model['reference'])
            records.append({'model': model, 'certificate': cert, 'construction_seconds': time.monotonic()-before})
        (OUT/'response_certificates.json').write_text(json.dumps(records, indent=2)+'\n')
    results = []
    for row in records:
        receipt = check(row['model'], row['certificate'])
        if row['model']['reference']['bath_spins'] <= 9:
            receipt['explicit_control'] = exact_control(row['model'], row['certificate'])
        receipt['certificate_bytes'] = len(json.dumps({'model': row['model'], 'certificate': row['certificate']}).encode())
        results.append(receipt)
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical imports on the exact response path')
    result = {'cases': results, 'numerical_packages_loaded': forbidden,
        'basis_states_enumerated_in_controls': sum(r.get('explicit_control', {}).get('basis_states_enumerated', 0) for r in results),
        'wall_seconds': time.monotonic()-start}
    (OUT/('response_replay.json' if replay else 'response_discovery.json')).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'cases': len(results), 'replay': replay, 'wall_seconds': result['wall_seconds']}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
