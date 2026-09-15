"""Density-square tangents with an orbital-sized exact fermionic filling bound."""
from fractions import Fraction as F
import time

from research.compact_response_20260913 import program
from research.molecular_collective_20260913.core import extract, tail_replay, digest
from research.positive_response_20260913.coercivity import spatial_inputs
from research.positive_response_20260913.molecular_diagnostic import spatial_one_body
from research.certificate_scaling.commutator_dual_witness import psd


def ingredients(data, tail, orbital, shifts):
    s = data['modes']//2; n = data['particles']
    if data['modes'] not in (12, 16) or n != s or type(orbital) is not int or not 0 <= orbital < s:
        raise ValueError('Bounded half-filled molecular sector required')
    p = extract(data, tail['center_number']); constant, t = spatial_one_body(p)
    matrices, weights, _ = spatial_inputs(data, tail); tr = tail_replay(data, tail)
    if len(shifts) != len(weights) or any(not isinstance(x, F) for x in shifts):
        raise ValueError('One exact tangent per density factor required')
    rem = [i for i in range(s) if i != orbital]; j = orbital
    V = [[sum(w*L[i][j]*L[k][j] for w, L in zip(weights, matrices)) for k in rem] for i in rem]
    leakage = sum(w*sum(L[i][j]**2 for i in rem) for w, L in zip(weights, matrices))
    effective = [[t[i][k]-V[a][b]/2+sum(w*x*L[i][k] for w, x, L in zip(weights, shifts, matrices))
                  for b, k in enumerate(rem)] for a, i in enumerate(rem)]
    base = constant+2*t[j][j]+leakage+F(tr['lower_operator_shift_Ha'])
    base += sum(w*(2*L[j][j]*x-x*x/2) for w, x, L in zip(weights, shifts, matrices))
    return base, effective, n-2, matrices, weights, t


def check(data, tail, cert):
    start = time.monotonic()
    if cert.get('kind') != 'spatial_tangent_sector_gap_v1' or cert['fixture_sha256'] != digest(data) or cert['tail_sha256'] != digest(tail):
        raise ValueError('Spatial gap binding failed')
    values = cert['tangents']
    if any(type(x) is not str for x in values) or type(cert['chemical_potential']) is not str:
        raise ValueError('Rational spatial certificate required')
    base, S, n, _, _, _ = ingredients(data, tail, cert['orbital'], list(map(F, values)))
    d = len(S); tau = F(cert['chemical_potential']); raw = cert['negative_part_majorant']
    if len(raw) != d or any(len(row) != d or any(type(x) is not str for x in row) for row in raw):
        raise ValueError('Wrong orbital majorant dimensions')
    Y = [list(map(F, row)) for row in raw]
    psd(Y); psd([[S[i][j]+Y[i][j]-(tau if i == j else 0) for j in range(d)] for i in range(d)])
    # S >= tau I-Y and spin occupancy <=2 imply dGamma(S)>=tau*N-2 tr(Y).
    lower = base+n*tau-2*sum(Y[i][i] for i in range(d))
    if type(cert['lower_H_QQ_Ha']) is not str or F(cert['lower_H_QQ_Ha']) > lower:
        raise ValueError('Spatial sector endpoint overstated')
    s=data['modes']//2;coefficient_dimension=s*(s+1)//2
    return {'lower_H_QQ_Ha': cert['lower_H_QQ_Ha'], 'proved_lower_H_QQ_Ha': str(lower),
        'largest_matrix_dimension': coefficient_dimension, 'new_orbital_matrix_dimension': d,
        'tail_coefficient_matrix_dimension':coefficient_dimension,
        'new_orbital_PSD_checks': 2, 'density_tangents': len(values),
        'many_body_states_enumerated': 0, 'many_body_matrix_entries': 0,
        'replay_seconds': time.monotonic()-start}


def propose(data, tail, orbital):
    import numpy as np
    from scipy.optimize import minimize
    start = time.monotonic(); matrices, weights, _ = spatial_inputs(data, tail)
    zeros = [F(0)]*len(weights); base, S, n, _, _, _ = ingredients(data, tail, orbital, zeros)
    rem = [i for i in range(data['modes']//2) if i != orbital]
    A = np.array([[[float(L[i][j]) for j in rem] for i in rem] for L in matrices])
    w = np.array(list(map(float, weights))); ell = np.array([float(L[orbital][orbital]) for L in matrices]); s = np.array(S, dtype=float)
    def objective(x):
        e, U = np.linalg.eigh(s+np.einsum('k,kij->ij', w*x, A)); occ = U[:, :n//2]
        value = float(base)+np.sum(w*(2*ell*x-x*x/2))+2*np.sum(e[:n//2])
        gradient = w*(2*ell-x+2*np.einsum('ia,kij,ja->k', occ, A, occ))
        return -value, -gradient
    result = minimize(objective, np.zeros(len(w)), jac=True, method='BFGS', options={'maxiter': 250, 'gtol': 1e-9})
    shifts = [F(round(float(x)*10**8), 10**8) for x in result.x]
    base, S, n, _, _, _ = ingredients(data, tail, orbital, shifts)
    e, U = np.linalg.eigh(np.array(S, dtype=float)); k = n//2
    tau = F(round(float((e[k-1]+e[k])/2)*10**8), 10**8)
    Yf = (U*np.maximum(float(tau)-e, 0.))@U.T+1e-8*np.eye(len(S))
    Y = [[F(round(float(x)*10**10), 10**10) for x in row] for row in Yf]
    lower = base+n*tau-2*sum(Y[i][i] for i in range(len(S)))
    cert = {'kind': 'spatial_tangent_sector_gap_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'orbital': orbital, 'tangents': list(map(str, shifts)), 'chemical_potential': str(tau),
        'negative_part_majorant': [[str(x) for x in row] for row in Y],
        'lower_H_QQ_Ha': str(program.floor_grid(lower, 10**10))}
    receipt = check(data, tail, cert)
    return cert, {'construction_seconds': time.monotonic()-start, 'iterations': int(result.nit),
        'optimizer_success': bool(result.success), 'optimizer_message': str(result.message),
        'acceptance': receipt, 'many_body_states_enumerated': 0}
