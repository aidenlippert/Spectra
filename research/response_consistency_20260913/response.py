"""Scalar exact certificates for the homogeneous central-spin Hamiltonian only."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/response_consistency_20260913'


def parameters(model):
    if set(model) != {'kind', 'bath_spins', 'excitations', 'epsilon', 'omega', 'chi', 'g'}:
        raise ValueError('Only the specified homogeneous Hamiltonian is supported')
    if model['kind'] != 'homogeneous_central_spin_v1':
        raise ValueError('Unknown response model')
    M = model['bath_spins']; q = model['excitations']
    if type(M) is not int or type(q) is not int or M < 1 or not 1 <= q <= M:
        raise ValueError('Nonempty retained and eliminated sectors required')
    if any(type(model[k]) is not str for k in ('epsilon', 'omega', 'chi', 'g')):
        raise ValueError('Exact rational model parameters required')
    epsilon, omega, chi, g = (F(model[k]) for k in ('epsilon', 'omega', 'chi', 'g'))
    if not g:
        raise ValueError('This coupled response control requires nonzero g')
    return M, q, omega*q, epsilon+(omega+chi)*(q-1), g, q*(M-q+1)


def check(model, cert):
    start = time.monotonic(); M, q, A0, D0, g, c = parameters(model)
    if set(cert) != {'lower', 'trial_t'} or any(type(x) is not str for x in cert.values()):
        raise ValueError('Exact lower and trial response coefficients required')
    b = F(cert['lower']); t = F(cert['trial_t']); d = D0-b
    if d <= 0:
        raise ValueError('Strict eliminated-sector gap required')
    slack = (A0-b)*d-g*g*c
    if slack < 0:
        raise ValueError('Collective Schur remainder is not positive')
    U = (A0+t*t*c*D0+2*g*t*c)/(1+t*t*c)
    if U < b:
        raise AssertionError('Analytic trial conflicts with the lower certificate')
    return {'lower': str(b), 'upper': str(U), 'width': str(U-b), 'width_float': float(U-b),
        'eliminated_gap': str(d), 'response_coefficient': str(g/d), 'residual_E': '0',
        'scalar_determinant_slack': str(slack), 'collective_norm_squared': c,
        'many_body_states_enumerated': 0, 'many_body_inverse_entries': 0,
        'replay_seconds': time.monotonic()-start,
        'scope': 'Exact ground interval only for the formula-specified homogeneous central-spin model in the fixed-q sector.'}


def propose(model, steps=48):
    _, _, A0, D0, g, c = parameters(model)
    # Since c>=1, |g|c bounds |g|sqrt(c). This brackets the smaller root.
    lo = min(A0, D0)-abs(g)*c-1; hi = min(A0, D0)
    for _ in range(steps):
        mid = (lo+hi)/2
        if (A0-mid)*(D0-mid) >= g*g*c:
            lo = mid
        else:
            hi = mid
    return {'lower': str(lo), 'trial_t': str(-g/(D0-lo))}


def cases():
    return [{'kind': 'homogeneous_central_spin_v1', 'bath_spins': M, 'excitations': M//2,
             'epsilon': epsilon, 'omega': '0', 'chi': '0', 'g': str(F(1, M))}
            for M, epsilon in ((2, '3/10'), (4, '3/10'), (6, '3/10'), (8, '3/10'),
                               (16, '-2/5'), (64, '0'), (256, '7/10'), (1024, '-1/100'),
                               (1000000, '3/10'))]


def run():
    start = time.monotonic(); records = []
    for model in cases():
        before = time.monotonic(); cert = propose(model); discovery = time.monotonic()-before
        accepted = check(model, cert)
        records.append({'model': model, 'certificate': cert, 'accepted': accepted,
            'bisection_seconds': discovery, 'certificate_bytes': len(json.dumps({'model': model, 'certificate': cert}).encode()),
            'maximum_rational_bits': max(max(F(x).numerator.bit_length(), F(x).denominator.bit_length()) for x in cert.values()),
            'physical_sector_dimension_formula': f"binomial({model['bath_spins']},{model['excitations']})+binomial({model['bath_spins']},{model['excitations']-1})"})
    forbidden = [n for n in ('numpy', 'scipy', 'cvxpy', 'pyscf') if n in sys.modules]
    if forbidden:
        raise AssertionError('Numerical imports on the exact response path')
    result = {'cases': records, 'numerical_packages_loaded': forbidden, 'wall_seconds': time.monotonic()-start,
        'scope': 'Known solvable family and proof-rule control; no generic molecular scaling claim.'}
    (OUT/'response_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
