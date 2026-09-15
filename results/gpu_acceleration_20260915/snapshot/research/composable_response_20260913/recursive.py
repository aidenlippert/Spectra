"""Second-level operator envelope and an exact obstruction to its scalar data."""
from fractions import Fraction as F
import argparse
import json
import sys
import time

from research.compact_response_20260913 import program
from research.molecular_collective_20260913.core import digest
from research.composable_response_20260913 import spatial_gap, coupling

OUT = program.ROOT/'results/composable_response_20260913'


def scalar_polynomial(d, response):
    delta = F(response['delta_Ha']); M = F(response['M_Ha']); c = (M+delta)/2; a = (M-delta)/2
    residual = program.chebyshev(response['order'], (c-d)/a)/program.chebyshev(response['order'], c/a)
    return (1-residual)/d


def counterexample(alpha, delta, M, nu, response):
    a = program.ceil_grid(max(alpha, F(1, 100000)), 10**5)
    d = program.ceil_grid(delta, 10**5)
    t = coupling.sqrt_up(a*d, 10**5)+F(1, 1000)
    if d > M or t*t > nu or t > F(response['coupling_norm_Ha']):
        raise ValueError('The bounded scalar counterexample does not fit this envelope')
    p = scalar_polynomial(d, response); K = a-t*t*(2*p-d*p*p)
    if K >= 0:
        raise ValueError('Scalar information obstruction not established')
    return {'a_Ha': str(a), 'd_Ha': str(d), 'coupling_Ha': str(t),
        'retained_K_Ha': str(K), 'witness_coordinates': ['1', str(-t*p)],
        'scope': 'A two-dimensional operator with the same scalar lower/upper/norm information makes the actual prescribed polynomial retained operator negative. This refutes sufficiency of those scalar data, not positivity of the actual molecular sector.'}


def check_counterexample(cert, alpha, delta, M, nu, response):
    a = F(cert['a_Ha']); d = F(cert['d_Ha']); t = F(cert['coupling_Ha'])
    if a < alpha or not delta <= d <= M or t*t > nu or abs(t) > F(response['coupling_norm_Ha']):
        raise ValueError('Counterexample violates the supplied information')
    p = scalar_polynomial(d, response); K = a-2*t*t*p+d*t*t*p*p
    if K >= 0 or str(K) != cert['retained_K_Ha'] or cert['witness_coordinates'] != ['1', str(-t*p)]:
        raise ValueError('Counterexample negative expectation did not reproduce')
    return {'exact_retained_K_Ha': str(K), 'retained_K_float_Ha': float(K), 'toy_dimension': 2}


def check(data, tail, response, cert):
    start = time.monotonic(); r = program.check(data, tail, response)
    if cert.get('kind') != 'second_elimination_envelope_v1' or cert['fixture_sha256'] != digest(data) or cert['response_sha256'] != digest(response):
        raise ValueError('Recursive-envelope input binding failed')
    s = data['modes']//2
    if [c['orbital'] for c in cert['sector_gaps']] != [s-1, s-2]:
        raise ValueError('Wrong nested occupation sectors')
    gaps = [spatial_gap.check(data, tail, c) for c in cert['sector_gaps']]
    b = F(response['target_Ha']); delta = F(gaps[0]['lower_H_QQ_Ha'])-b; alpha = F(gaps[1]['lower_H_QQ_Ha'])-b
    if delta <= 0:
        raise ValueError('No positive first-sector endpoint')
    c = coupling.prove(data); nu = F(c['coupling_norm_squared_Ha2']); eta = F(r['exact_program_residual_penalty_Ha'])
    bound = alpha-nu/delta-eta
    if F(cert['second_gap_lower_Ha']) != bound:
        raise ValueError('Second-sector envelope did not reproduce')
    counter = check_counterexample(cert['scalar_information_obstruction'], alpha, delta, F(response['M_Ha']), nu, response) if bound <= 0 else None
    return {'first_gap_Ha': str(delta), 'first_gap_float_Ha': float(delta),
        'second_bare_gap_Ha': str(alpha), 'second_bare_gap_float_Ha': float(alpha),
        'induced_correction_upper_Ha': str(nu/delta), 'induced_correction_upper_float_Ha': float(nu/delta),
        'response_penalty_Ha': str(eta), 'second_gap_lower_Ha': str(bound), 'second_gap_lower_float_Ha': float(bound),
        'second_positive_gap_certified': bound > 0, 'scalar_information_obstruction': counter,
        'spatial_gap_checks': gaps, 'coupling': c, 'many_body_states_enumerated': 0,
        'many_body_matrix_entries': 0, 'H_actions_per_retained_operator_application_with_shared_initial_action': 2*response['order']+1,
        'replay_seconds': time.monotonic()-start,
        'scope': 'Constructed from CAR words and orbital matrices. A failed scalar gap estimate is inconclusive for the actual molecular second gap; the counterexample isolates the missing correlation information.'}


def propose(data, tail, response):
    start = time.monotonic(); program_receipt = program.check(data, tail, response); b = F(response['target_Ha'])
    gaps = []; discovery = []
    for j in (data['modes']//2-1, data['modes']//2-2):
        cert, stats = spatial_gap.propose(data, tail, j); gaps.append(cert); discovery.append(stats)
    c = coupling.prove(data); delta = F(gaps[0]['lower_H_QQ_Ha'])-b; alpha = F(gaps[1]['lower_H_QQ_Ha'])-b
    nu = F(c['coupling_norm_squared_Ha2']); eta = F(program_receipt['exact_program_residual_penalty_Ha']); bound = alpha-nu/delta-eta
    obstruction = counterexample(alpha, delta, F(response['M_Ha']), nu, response) if bound <= 0 else None
    cert = {'kind': 'second_elimination_envelope_v1', 'fixture_sha256': digest(data),
        'response_sha256': digest(response), 'sector_gaps': gaps, 'second_gap_lower_Ha': str(bound),
        'scalar_information_obstruction': obstruction}
    return cert, {'spatial_gap_discovery': discovery, 'coupling': c, 'construction_seconds': time.monotonic()-start}


def affine(*terms):
    return [sum((coefficient*vector[i] for coefficient, vector in terms), F(0)) for i in range(len(terms[0][1]))]


def apply_p(action_D, vector, response):
    delta = F(response['delta_Ha']); M = F(response['M_Ha']); c = (M+delta)/2; a = (M-delta)/2; z = c/a
    previous = [F(0)]*len(vector); current = affine((1/a, vector)); tp = F(1); tc = z
    for _ in range(1, response['order']):
        Z = affine((c/a, current), (-1/a, action_D(current)))
        previous, current = current, affine((F(2), Z), (F(-1), previous), (2*tc/a, vector))
        tp, tc = tc, 2*z*tc-tp
    return affine((1/tc, current))


def retained_action(action_shifted_H, project_P, project_Q, vector, response):
    """Composable oracle action for K, with first H action shared by A and C.

    The caller owns representation and action costs. This function supplies
    no positivity certificate and does not assume an oracle action is cheap.
    """
    Hv = action_shifted_H(project_P(vector)); W = project_Q(Hv)
    D = lambda x: project_Q(action_shifted_H(project_Q(x)))
    pW = apply_p(D, W, response); DpW = D(pW); pDpW = apply_p(D, DpW, response)
    FW = affine((F(2), pW), (F(-1), pDpW))
    return affine((F(1), project_P(Hv)), (F(-1), project_P(action_shifted_H(FW))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true'); replay = parser.parse_args().replay
    start = time.monotonic(); rows = []
    for name in ('h6', 'h8'):
        data, tail, reference = program.load_case(name); response = json.loads((program.OUT/f'{name}_program.json').read_text())
        path = OUT/f'{name}_recursive.json'; discovery = None
        if replay:
            cert = json.loads(path.read_text())
        else:
            cert, discovery = propose(data, tail, response); path.write_text(json.dumps(cert, indent=2)+'\n')
        receipt = check(data, tail, response, cert)
        rows.append({'case': name, 'discovery': discovery, 'receipt': receipt, 'certificate_bytes': path.stat().st_size})
    forbidden = [x for x in ('numpy', 'scipy', 'cvxpy', 'pyscf') if x in sys.modules]
    if replay and forbidden:
        raise AssertionError('Numerical imports during recursive replay')
    result = {'cases': rows, 'wall_seconds': time.monotonic()-start, 'numerical_packages_loaded': forbidden}
    (OUT/('recursive_replay.json' if replay else 'recursive_discovery.json')).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'cases': [(r['case'], r['receipt']['second_gap_lower_float_Ha']) for r in rows], 'wall_seconds': result['wall_seconds']}), flush=True)
