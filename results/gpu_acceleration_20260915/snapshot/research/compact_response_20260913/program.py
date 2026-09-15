"""A bounded Chebyshev response program and exact molecular endpoint proofs."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys
import time

from experiments.marginal_symbolic import add, mono, product, scale
from research.molecular_collective_20260913.core import digest, extract, tail_replay
from research.positive_response_20260913.molecular_diagnostic import spatial_one_body
from research.positive_response_20260913.coercivity import spatial_inputs, spatial_norm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/compact_response_20260913'
KIND = 'chebyshev_double_occupancy_response_v1'


def load_case(name):
    if name not in ('h6', 'h8'):
        raise ValueError('This experiment is capped at frozen H6 and H8')
    rank = {'h6': 10, 'h8': 14}[name]
    src = ROOT/'results/molecular_collective_20260913/campaign'/name
    data = json.loads((src/'fixture.json').read_text()); tail = json.loads((src/f'rank_{rank}/tail.json').read_text())
    upper_path = ROOT/'results/certificate_scaling/active_space_ladder_references_aligned'/name/'upper.json'
    return data, tail, json.loads(upper_path.read_text())


def floor_grid(value, den=10000):
    return F(value.numerator*den//value.denominator, den)


def ceil_grid(value, den=10000):
    return -floor_grid(-value, den)


def projector(m):
    return product(mono(((1, m-2), (0, m-2))), mono(((1, m-1), (0, m-1))))


def one_body_endpoint(matrix, n, upper=False):
    if n % 2 or not 0 < n < 2*len(matrix):
        raise ValueError('Even nontrivial remaining particle count required')
    rho = max(sum(abs(x) for j, x in enumerate(row) if i != j) for i, row in enumerate(matrix))
    diagonal = sorted(matrix[i][i] for i in range(len(matrix)))
    return 2*sum(diagonal[-n//2:] if upper else diagonal[:n//2])+(n*rho if upper else -n*rho)


def bounds(data, tail, target):
    start = time.monotonic(); m = data['modes']; n = data['particles']
    if m not in (12, 16) or n != m//2:
        raise ValueError('This experiment requires the half-filled H6 or H8 sector')
    p = extract(data, tail['center_number']); matrices, weights, _ = spatial_inputs(data, tail)
    constant, t = spatial_one_body(p); tail_receipt = tail_replay(data, tail)
    last = m//2-1; nr = n-2
    V = [[F(0) for _ in range(last)] for _ in range(last)]
    leakage_constant = F(0); interaction_upper = F(0)
    for matrix, weight in zip(matrices, weights):
        A = [row[:last] for row in matrix[:last]]; v = [matrix[i][last] for i in range(last)]; ell = matrix[last][last]
        vv = sum(x*x for x in v); leakage_constant += weight*vv
        interaction_upper += weight*((nr*spatial_norm(A)+2*abs(ell))**2/2+vv)
        for i in range(last):
            for j in range(last):
                V[i][j] += weight*v[i]*v[j]
    effective = [[t[i][j]-V[i][j]/2 for j in range(last)] for i in range(last)]
    core = [row[:last] for row in t[:last]]; base = constant+2*t[last][last]
    gamma = base+leakage_constant+one_body_endpoint(effective, nr)+F(tail_receipt['lower_operator_shift_Ha'])
    gamma_old = base+one_body_endpoint(core, nr)+F(tail_receipt['lower_operator_shift_Ha'])
    upper = base+one_body_endpoint(core, nr, True)+interaction_upper+F(tail_receipt['upper_operator_shift_Ha'])
    Q = projector(m); P = add(mono(()), scale(Q, -1)); coupling = product(product(Q, p['h']), P)
    g = sum(abs(c) for c in coupling.values())
    return {'delta': gamma-target, 'M': upper-target, 'coupling_norm': g,
        'gamma': gamma, 'upper_H_QQ': upper, 'old_delta': gamma_old-target,
        'leakage_constant': leakage_constant, 'tail': tail_receipt,
        'coupling_operator_terms': len(coupling), 'coupling_max_degree': max(map(len, coupling), default=0),
        'many_body_states_enumerated': 0, 'wall_seconds': time.monotonic()-start}


def chebyshev(k, z):
    if type(k) is not int or not 0 <= k <= 256:
        raise ValueError('Chebyshev order must lie in 0..256')
    previous = F(1)
    if k == 0:
        return previous
    current = z
    for _ in range(1, k):
        previous, current = current, 2*z*current-previous
    return current


def propose(data, tail, target, budget=F(1, 1000000)):
    start = time.monotonic(); interval = bounds(data, tail, target)
    delta = floor_grid(interval['delta']); M = ceil_grid(interval['M']); g = ceil_grid(interval['coupling_norm'])
    if not 0 < delta < M:
        raise ValueError('Independent strictly positive nonscalar spectral interval required')
    z = (M+delta)/(M-delta)
    for k in range(1, 257):
        Tk = chebyshev(k, z)
        if g*g <= budget*delta*Tk*Tk:
            break
    else:
        raise ValueError('Response-order budget exhausted')
    cert = {'kind': KIND, 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'doubly_occupied_spatial_orbital': data['modes']//2-1, 'target_Ha': str(target),
        'delta_Ha': str(delta), 'M_Ha': str(M), 'coupling_norm_Ha': str(g),
        'order': k, 'residual_penalty_budget_Ha': str(budget)}
    return cert, {'construction_seconds': time.monotonic()-start, 'bounds_seconds': interval['wall_seconds'],
        'many_body_states_enumerated': 0, 'many_body_inverse_entries': 0}


def check(data, tail, cert):
    start = time.monotonic()
    keys = {'kind', 'fixture_sha256', 'tail_sha256', 'doubly_occupied_spatial_orbital', 'target_Ha',
            'delta_Ha', 'M_Ha', 'coupling_norm_Ha', 'order', 'residual_penalty_budget_Ha'}
    if set(cert) != keys or cert['kind'] != KIND:
        raise ValueError('Unknown response program')
    if cert['fixture_sha256'] != digest(data) or cert['tail_sha256'] != digest(tail):
        raise ValueError('Response input binding failed')
    orbital = cert['doubly_occupied_spatial_orbital']
    if type(orbital) is not int or orbital != data['modes']//2-1:
        raise ValueError('Wrong eliminated double-occupancy projector')
    names = ('target_Ha', 'delta_Ha', 'M_Ha', 'coupling_norm_Ha', 'residual_penalty_budget_Ha')
    if any(type(cert[key]) is not str for key in names):
        raise ValueError('Exact rational response parameters required')
    target, delta, M, g, budget = (F(cert[key]) for key in names); k = cert['order']
    if type(k) is not int or not 1 <= k <= 256 or not 0 < delta < M or min(g, budget) <= 0:
        raise ValueError('Invalid response order, spectral interval or budget')
    interval = bounds(data, tail, target)
    if delta > interval['delta'] or M < interval['M'] or g < interval['coupling_norm']:
        raise ValueError('Claimed spectral or coupling bound was not independently proved')
    z = (M+delta)/(M-delta); Tk = chebyshev(k, z); eta = g*g/(delta*Tk*Tk)
    if Tk <= 1 or eta > budget:
        raise ValueError('Response residual penalty exceeds the claimed budget')
    return {'response_program_sha256': digest(cert), 'target_Ha': str(target),
        'order': k, 'polynomial_degree': k-1, 'D_actions_per_right_hand_side': k-1,
        'exact_program_residual_factor': str(1/Tk), 'residual_factor_float': float(1/Tk),
        'exact_program_residual_penalty_Ha': str(eta), 'residual_penalty_float_Ha': float(eta),
        'delta_Ha': str(delta), 'delta_float_Ha': float(delta), 'M_Ha': str(M), 'M_float_Ha': float(M),
        'coupling_norm_Ha': str(g), 'old_delta_float_Ha': float(interval['old_delta']),
        'unrounded_delta_float_Ha': float(interval['delta']),
        'coupling_operator_terms': interval['coupling_operator_terms'],
        'coupling_max_degree': interval['coupling_max_degree'],
        'scalar_recurrence_rational_bits': max(Tk.numerator.bit_length(), Tk.denominator.bit_length()),
        'many_body_states_enumerated': 0, 'many_body_inverse_entries': 0,
        'retained_remainder_certified_by_this_rule': False,
        'replay_seconds': time.monotonic()-start,
        'scope': 'Exact operator-program error certificate on every retained vector. K positivity and any response-action cost are separate.'}


def apply_response(action_D, W, cert, order=None):
    """Numerical proposer only: array W may contain one or many right-hand sides."""
    k = cert['order'] if order is None else order
    delta = float(F(cert['delta_Ha'])); M = float(F(cert['M_Ha'])); c = (M+delta)/2; a = (M-delta)/2; z = c/a
    q_previous = 0*W; q_current = W/a; T_previous = 1.; T_current = z
    for _ in range(1, k):
        Zq = (c*q_current-action_D(q_current))/a
        q_previous, q_current = q_current, 2*Zq-q_previous+(2/a)*T_current*W
        T_previous, T_current = T_current, 2*z*T_current-T_previous
    return q_current/T_current


def run(replay=False):
    start = time.monotonic(); OUT.mkdir(exist_ok=True); rows = []
    for name in ('h6', 'h8'):
        data, tail, reference = load_case(name); path = OUT/f'{name}_program.json'
        if replay:
            cert = json.loads(path.read_text()); discovery = None
        else:
            cert, discovery = propose(data, tail, F(reference['upper'])-F(1, 1000))
            path.write_text(json.dumps(cert, indent=2)+'\n')
        receipt = check(data, tail, cert); receipt['molecule'] = name; receipt['certificate_bytes'] = path.stat().st_size
        if discovery is not None:
            receipt['discovery'] = discovery
        rows.append(receipt)
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical imports on the exact response-program path')
    result = {'cases': rows, 'numerical_packages_loaded': forbidden, 'wall_seconds': time.monotonic()-start}
    (OUT/('program_replay.json' if replay else 'program_discovery.json')).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'cases': [(r['molecule'], r['order'], r['delta_float_Ha'], r['residual_penalty_float_Ha']) for r in rows],
        'wall_seconds': result['wall_seconds']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
