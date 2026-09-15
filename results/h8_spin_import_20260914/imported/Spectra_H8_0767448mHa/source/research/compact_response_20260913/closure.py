"""Explicit H6 retained-closure control; every sector and matrix is charged."""
from fractions import Fraction as F
from itertools import combinations
from math import lcm
import argparse
import json
import sys
import time

from experiments.marginal_symbolic import decode
from research.certificate_scaling.streaming_reference_upper import compile_term, upper
from research.molecular_collective_20260913.core import digest
from research.compact_response_20260913 import program

OUT = program.OUT


def blocks(data):
    start = time.monotonic()
    if data['modes'] != 12 or data['particles'] != 6:
        raise ValueError('Explicit retained closure is capped at H6')
    h = decode(data['hamiltonian'], 12, 4); den = lcm(*(c.denominator for c in h.values()))
    parity = sum(3 << (2*i) for i in range(0, 6, 2)); up_mask = sum(1 << (2*i) for i in range(6))
    terms = []
    for word, coefficient in h.items():
        flip = 0
        for _, mode in word:
            flip ^= 1 << mode
        if (flip & parity).bit_count() % 2 or sum((2*c-1) for c, mode in word if mode % 2 == 0):
            raise ValueError('Hamiltonian does not preserve the stated complete block partition')
        term = compile_term(word, int(coefficient*den))
        if term is not None:
            terms.append(term)
    groups = {}
    for occupied in combinations(range(12), 6):
        state = sum(1 << i for i in occupied); key = ((state & up_mask).bit_count(), (state & parity).bit_count() % 2)
        groups.setdefault(key, []).append(state)
    out = []
    for key, states in sorted(groups.items()):
        positions = {state: i for i, state in enumerate(states)}; H = [[0]*len(states) for _ in states]
        for j, state in enumerate(states):
            for required, occupied, flip, parity_term, coefficient in terms:
                if state & required == occupied:
                    target = state ^ flip
                    if target not in positions:
                        raise ValueError('Hamiltonian action left a claimed conserved block')
                    H[positions[target]][j] += coefficient*(-1 if (state & parity_term).bit_count() % 2 else 1)
        if any(H[i][j] != H[j][i] for i in range(len(H)) for j in range(i)):
            raise ValueError('Exact block matrix is not symmetric')
        p = [i for i, state in enumerate(states) if state & (3 << 10) != 3 << 10]
        q = [i for i, state in enumerate(states) if state & (3 << 10) == 3 << 10]
        out.append({'key': list(key), 'states': states, 'H': H, 'p': p, 'q': q})
    return out, den, {'physical_sector_basis_labels_generated': sum(len(g['states']) for g in out),
        'Hamiltonian_action_sources': sum(len(g['states']) for g in out),
        'word_state_checks': sum(len(g['states'])*len(terms) for g in out),
        'retained_dimension': sum(len(g['p']) for g in out), 'eliminated_dimension': sum(len(g['q']) for g in out),
        'largest_full_conserved_block': max(len(g['states']) for g in out),
        'largest_retained_block': max(len(g['p']) for g in out),
        'largest_eliminated_block': max(len(g['q']) for g in out),
        'full_block_matrix_entries': sum(len(g['states'])**2 for g in out),
        'response_matrix_entries': sum(len(g['p'])*len(g['q']) for g in out),
        'retained_matrix_entries': sum(len(g['p'])**2 for g in out),
        'block_partition_parity_mask': parity, 'preparation_seconds': time.monotonic()-start}


def shifted(block, hden, target):
    den = lcm(hden, target.denominator); factor = den//hden; shift = int(target*den)
    H = block['H']; p = block['p']; q = block['q']
    A = [[H[i][j]*factor-(shift if i == j else 0) for j in p] for i in p]
    D = [[H[i][j]*factor-(shift if i == j else 0) for j in q] for i in q]
    Bt = [[H[i][j]*factor for j in p] for i in q]
    return A, Bt, D, den


def residual_and_K(A, Bt, D, den, X, xden, delta):
    n = len(A); q = len(D)
    if len(X) != q or any(len(row) != n or any(type(x) is not int for x in row) for row in X):
        raise ValueError('Wrong exact response dimensions')
    rows = [[(k, d) for k, d in enumerate(row) if d] for row in D]
    Y = [[sum(d*X[k][j] for k, d in row) for j in range(n)] for row in rows]
    e2 = sum((Bt[i][j]*xden-Y[i][j])**2 for i in range(q) for j in range(n))
    eta = F(e2, den*den*xden*xden)/delta
    Cx = list(zip(*X)) if q else [()]*n; Cb = list(zip(*Bt)) if q else [()]*n; Cy = list(zip(*Y)) if q else [()]*n
    K = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            value = A[i][j]*xden*xden-xden*(sum(x*b for x, b in zip(Cx[i], Cb[j]))+sum(b*x for b, x in zip(Cb[i], Cx[j])))
            value += sum(x*y for x, y in zip(Cx[i], Cy[j]))
            K[i][j] = K[j][i] = value
    return K, den*xden*xden, eta


def check_factor(K, kden, eta, factor, fden):
    n = len(K)
    if len(factor) != n or any(len(row) != i+1 or any(type(x) is not int for x in row) for i, row in enumerate(factor)):
        raise ValueError('A complete exact triangular retained factor is required')
    den = lcm(kden, eta.denominator, fden*fden); km = den//kden; fm = den//(fden*fden); e = int(eta*den)
    rows = [0]*n
    for i in range(n):
        for j in range(i+1):
            residual = K[i][j]*km-sum(x*y for x, y in zip(factor[i], factor[j]))*fm-(e if i == j else 0)
            if i == j:
                rows[i] += residual
            else:
                rows[i] -= abs(residual); rows[j] -= abs(residual)
    if min(rows) < 0:
        raise ValueError('Retained remainder positivity did not pass exact residual bounds')
    return F(min(rows), den)


def check_scalar_obstruction(groups, hden, response, cert):
    if cert['response_program_sha256'] != digest(response):
        raise ValueError('Scalar obstruction binding failed')
    block = next((g for g in groups if g['key'] == cert['block_key']), None)
    if block is None:
        raise ValueError('Unknown scalar obstruction block')
    A, Bt, _, den = shifted(block, hden, F(response['target_Ha'])); v = cert['vector']
    if len(v) != len(A) or any(type(x) is not int for x in v) or not any(v):
        raise ValueError('Invalid scalar obstruction vector')
    norm = sum(x*x for x in v)
    a = F(sum(v[i]*A[i][j]*v[j] for i in range(len(v)) for j in range(len(v))), den)
    b2 = F(sum(sum(c*x for c, x in zip(row, v))**2 for row in Bt), den*den)
    value = (a-b2/F(response['delta_Ha']))/norm
    if value >= 0 or str(value) != cert['exact_expectation_Ha']:
        raise ValueError('Scalar-Schur failure did not reproduce')
    return {'exact_expectation_Ha': str(value), 'expectation_float_Ha': float(value),
        'scope': 'Refutes the scalar response B B*/delta at this target and the same independently certified delta.'}


def propose(data, tail, response):
    import numpy as np
    start = time.monotonic(); program.check(data, tail, response)
    groups, hden, cost = blocks(data); target = F(response['target_Ha']); delta = F(response['delta_Ha'])
    xden = 10**10; fden = 10**11; margin = F(1, 10000000)
    rows = []; orders = sorted({1, 2, 4, 8, 16, response['order']}); curve = {k: float('inf') for k in orders}; worst_scalar = None
    recurrence_actions = 0; last_X_actions = 0
    for group in groups:
        A, Bt, D, den = shifted(group, hden, target); n = len(A); q = len(D)
        af = np.array(A, dtype=float)/float(den); bf = np.array(Bt, dtype=float).reshape((q, n))/float(den)
        df = np.array(D, dtype=float).reshape((q, q))/float(den)
        scalar = af-bf.T@bf/float(delta); values, vectors = np.linalg.eigh(scalar)
        if worst_scalar is None or values[0] < worst_scalar[0]:
            worst_scalar = (float(values[0]), group, np.rint(vectors[:, 0]*10**9).astype(np.int64).tolist())
        for order in orders:
            Xf = program.apply_response(lambda x: df@x, bf, response, order) if q else np.zeros((0, n))
            recurrence_actions += (order-1)*n if q else 0
            E = bf-df@Xf; Kf = af-bf.T@Xf-Xf.T@bf+Xf.T@df@Xf
            eta_float = float(np.sum(E*E))/float(delta)
            curve[order] = min(curve[order], float(np.linalg.eigvalsh(Kf)[0])-eta_float)
            if order == response['order']:
                selected = Xf
        last_X_actions += (response['order']-1)*n if q else 0
        X = np.rint(selected*xden).astype(np.int64).tolist()
        K, kden, eta_exact = residual_and_K(A, Bt, D, den, X, xden, delta)
        eta = program.ceil_grid(eta_exact, 10**12)
        kf = np.array(K, dtype=float)/float(kden)
        L = np.linalg.cholesky(kf-(float(eta)+float(margin))*np.eye(n))
        factor = [[int(round(float(L[i, j])*fden)) for j in range(i+1)] for i in range(n)]
        checked = check_factor(K, kden, eta, factor, fden)
        rows.append({'key': group['key'], 'response_columns': X, 'retained_factor': factor,
            'residual_penalty_Ha': str(eta), 'checked_retained_margin_Ha': str(checked)})
    _, group, v = worst_scalar; A, Bt, _, den = shifted(group, hden, target)
    norm = sum(x*x for x in v); a = F(sum(v[i]*A[i][j]*v[j] for i in range(len(v)) for j in range(len(v))), den)
    b2 = F(sum(sum(c*x for c, x in zip(row, v))**2 for row in Bt), den*den)
    obstruction = {'response_program_sha256': digest(response), 'block_key': group['key'], 'vector': v,
        'exact_expectation_Ha': str((a-b2/delta)/norm)}
    check_scalar_obstruction(groups, hden, response, obstruction)
    cert = {'kind': 'expanded_h6_response_closure_v1', 'fixture_sha256': digest(data),
        'response_program_sha256': digest(response), 'response_denominator': xden,
        'retained_factor_denominator': fden, 'blocks': rows}
    stats = cost | {'construction_seconds': time.monotonic()-start,
        'all_degree_curve_D_actions_on_individual_columns': recurrence_actions,
        'selected_response_D_actions_on_individual_columns': last_X_actions,
        'numerical_degree_curve': [{'order': k, 'minimum_K_minus_Frobenius_residual_float_Ha': curve[k]} for k in orders],
        'scope': 'Explicit finite H6 retained-closure control. Response columns and retained factors are expanded and charged. No matrix inverse is used.'}
    return cert, obstruction, stats


def check(data, tail, response, cert, reference):
    start = time.monotonic(); response_receipt = program.check(data, tail, response)
    if cert.get('kind') != 'expanded_h6_response_closure_v1' or cert['fixture_sha256'] != digest(data) or cert['response_program_sha256'] != digest(response):
        raise ValueError('Full-closure binding failed')
    xden = cert['response_denominator']; fden = cert['retained_factor_denominator']
    if any(type(d) is not int or not 0 < d <= 10**14 for d in (xden, fden)):
        raise ValueError('Bounded positive exact denominators required')
    groups, hden, cost = blocks(data)
    if len(cert['blocks']) != len(groups) or [r['key'] for r in cert['blocks']] != [g['key'] for g in groups]:
        raise ValueError('Every conserved block must be covered exactly once')
    receipts = []; target = F(response['target_Ha']); delta = F(response['delta_Ha'])
    for group, row in zip(groups, cert['blocks']):
        A, Bt, D, den = shifted(group, hden, target)
        K, kden, actual_eta = residual_and_K(A, Bt, D, den, row['response_columns'], xden, delta)
        if type(row['residual_penalty_Ha']) is not str:
            raise ValueError('Rational residual penalty required')
        eta = F(row['residual_penalty_Ha'])
        if eta < actual_eta:
            raise ValueError('Rounded response residual penalty was understated')
        margin = check_factor(K, kden, eta, row['retained_factor'], fden)
        if str(margin) != row['checked_retained_margin_Ha']:
            raise ValueError('Retained margin did not reproduce')
        receipts.append({'key': group['key'], 'retained_dimension': len(group['p']), 'eliminated_dimension': len(group['q']),
            'actual_residual_penalty_Ha': str(actual_eta), 'accepted_residual_penalty_Ha': str(eta),
            'retained_margin_Ha': str(margin), 'retained_margin_float_Ha': float(margin)})
    U, upper_receipt = upper(data, reference['independent_upper'])
    if U != F(reference['upper']) or U < target:
        raise ValueError('Frozen upper did not reproduce or conflicts with the lower')
    return {'lower_Ha': str(target), 'upper_Ha': str(U), 'width_Ha': str(U-target), 'width_mHa': float(1000*(U-target)),
        'blocks': receipts, 'response_program_replay': response_receipt, 'upper_replay': upper_receipt,
        'explicit_enumeration_and_matrices': cost, 'matrix_inverse_entries': 0,
        'replay_seconds': time.monotonic()-start,
        'scope': 'A complete lower bound using explicitly expanded retained closure. The compact program alone does not prove this endpoint.'}


def run(replay=False):
    data, tail, reference = program.load_case('h6'); response = json.loads((OUT/'h6_program.json').read_text())
    if replay:
        start = time.monotonic(); cert = json.loads((OUT/'expanded_closure_certificate.json').read_text())
        result = check(data, tail, response, cert, reference)
        groups, hden, _ = blocks(data); obstruction = json.loads((OUT/'scalar_obstruction.json').read_text())
        result['scalar_obstruction'] = check_scalar_obstruction(groups, hden, response, obstruction)
        result['scalar_obstruction_extra_basis_labels'] = sum(len(g['states']) for g in groups)
        forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
        if forbidden:
            raise AssertionError('Numerical imports during exact expanded-closure replay')
        result['numerical_packages_loaded'] = forbidden; result['complete_replay_seconds'] = time.monotonic()-start
        (OUT/'closure_replay.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({'width_mHa': result['width_mHa'], 'replay_seconds': result['complete_replay_seconds']}), flush=True)
    else:
        cert, obstruction, stats = propose(data, tail, response)
        raw = json.dumps(cert, separators=(',', ':'))+'\n'; (OUT/'expanded_closure_certificate.json').write_text(raw)
        stats['expanded_certificate_bytes'] = len(raw.encode())
        (OUT/'scalar_obstruction.json').write_text(json.dumps(obstruction, indent=2)+'\n')
        (OUT/'closure_discovery.json').write_text(json.dumps(stats, indent=2)+'\n')
        print(json.dumps({'complete': True, 'construction_seconds': stats['construction_seconds'],
            'expanded_certificate_bytes': stats['expanded_certificate_bytes']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
