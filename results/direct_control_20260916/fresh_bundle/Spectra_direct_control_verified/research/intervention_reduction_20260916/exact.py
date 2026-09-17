"""Exact finite-sector Galerkin reduction with a uniform control envelope.

The only state labels generated are the supplied configurations and their CAR
neighbors. This can still exhaust a small sector; receipts report the counts.
No numerical dependency is imported in this accepting module.
"""
from fractions import Fraction as F
from itertools import product as cartesian_product
from math import comb, isqrt, lcm
import hashlib
import json
import time

from experiments.marginal_symbolic import decode, hermitian


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def rational(x):
    if not isinstance(x, str):
        raise ValueError('Rational inputs must be strings')
    return F(x)


def ceil_grid(x, scale=10**12):
    return F(-((-x.numerator * scale) // x.denominator), scale)


def sqrt_up(x, scale=10**12):
    if x < 0:
        raise ValueError('Negative square norm')
    n = isqrt(x.numerator * scale * scale // x.denominator)
    if F(n*n, scale*scale) < x:
        n += 1
    return F(n, scale)


def inverse_spd(a):
    """Exact positive LDL gate followed by exact inversion, no pivot tolerance."""
    n = len(a)
    if not n or any(len(row) != n for row in a):
        raise ValueError('Square metric required')
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Symmetric metric required')
    work = [list(map(F, row)) for row in a]
    for k in range(n):
        pivot = work[k][k]
        if pivot <= 0:
            raise ValueError('Metric is not positive definite')
        for i in range(k+1, n):
            for j in range(i, n):
                work[j][i] = work[i][j] = work[i][j] - work[i][k]*work[j][k]/pivot
    aug = [list(map(F, row)) + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        pivot = aug[k][k]
        aug[k] = [x/pivot for x in aug[k]]
        for i in range(n):
            if i != k:
                factor = aug[i][k]
                if factor:
                    aug[i] = [x-factor*y for x,y in zip(aug[i], aug[k])]
    return [row[n:] for row in aug]


def matmul(a, b):
    columns = list(zip(*b))
    return [[sum((x*y for x,y in zip(row, col)), F(0)) for col in columns] for row in a]


def apply_word(word, state):
    sign = 1
    for creation, mode in reversed(word):
        bit = 1 << mode
        if bool(state & bit) == bool(creation):
            return None
        if (state & (bit-1)).bit_count() % 2:
            sign = -sign
        state ^= bit
    return state, sign


def decode_model(data):
    m, n = data.get('modes'), data.get('particles')
    if type(m) is not int or type(n) is not int or not m > 0 or not 0 <= n <= m:
        raise ValueError('Invalid finite fermion sector')
    h = decode(data['hamiltonian'], m, 4)
    if not hermitian(h) or any(sum(2*c-1 for c,_ in w) for w in h):
        raise ValueError('Hermitian number-conserving Hamiltonian required')
    return h


def integer_action(poly, states):
    den = lcm(*(c.denominator for c in poly.values())) if poly else 1
    terms = [(w, int(c*den)) for w,c in poly.items()]
    columns = []
    for state in states:
        col = {}
        for word, coefficient in terms:
            result = apply_word(word, state)
            if result is not None:
                target, sign = result
                col[target] = col.get(target, 0) + sign*coefficient
        columns.append({s:v for s,v in col.items() if v})
    return den, columns


def weighted_action(columns, vectors):
    rank = len(vectors[0])
    result = {}
    for col, row in zip(columns, vectors):
        for state, coefficient in col.items():
            target = result.setdefault(state, [0]*rank)
            for j, value in enumerate(row):
                target[j] += coefficient*value
    return {s:row for s,row in result.items() if any(row)}


def exact_moments(data, proposal):
    if proposal.get('kind') != 'integer_control_subspace_v1' or proposal.get('fixture_sha256') != digest(data):
        raise ValueError('Proposal kind or fixture binding')
    h = decode_model(data)
    m, n = data['modes'], data['particles']
    states = proposal['configurations']
    if (not states or any(type(s) is not int or not 0 <= s < 1 << m or s.bit_count() != n for s in states)
            or len(set(states)) != len(states)):
        raise ValueError('Distinct valid fixed-N configurations required')
    vectors = proposal['vectors']
    den_v = proposal['denominator']
    if type(den_v) is not int or den_v <= 0 or len(vectors) != len(states) or not vectors or not vectors[0]:
        raise ValueError('Vector dimensions or denominator')
    rank = len(vectors[0])
    if any(len(row) != rank or any(type(v) is not int for v in row) for row in vectors):
        raise ValueError('Rectangular integer vectors required')
    controls = []
    amplitudes = []
    for item in proposal['controls']:
        p = decode(item['operator'], m, 4)
        if not hermitian(p) or any(sum(2*c-1 for c,_ in w) for w in p):
            raise ValueError('Hermitian number-conserving control required')
        amplitude = rational(item['amplitude_Ha'])
        if amplitude < 0:
            raise ValueError('Negative amplitude')
        controls.append(p)
        amplitudes.append(amplitude)
    if len(controls) > 4:
        raise ValueError('This implementation supports at most four controls')
    g = [[F(sum(row[i]*row[j] for row in vectors), den_v**2) for j in range(rank)] for i in range(rank)]
    inv_g = inverse_spd(g)
    raw_actions = [integer_action(poly, states) for poly in [h]+controls]
    actions = [weighted_action(cols, vectors) for _,cols in raw_actions]
    dens = [d for d,_ in raw_actions]
    ks = []
    for d, action in zip(dens, actions):
        k = [[F(sum(row[i]*action.get(s, [0]*rank)[j] for s,row in zip(states, vectors)), den_v**2*d)
              for j in range(rank)] for i in range(rank)]
        if any(k[i][j] != k[j][i] for i in range(rank) for j in range(rank)):
            raise AssertionError('Exact projected Hermiticity failed')
        ks.append(k)
    generators = [matmul(inv_g, k) for k in ks]
    tdiag = []
    for l, al in enumerate(actions):
        trow = []
        for q, aq in enumerate(actions):
            common = al.keys() & aq.keys()
            trow.append([F(sum(al[s][j]*aq[s][j] for s in common), den_v**2*dens[l]*dens[q]) for j in range(rank)])
        tdiag.append(trow)
    residuals = []
    for l in range(len(actions)):
        row = []
        for q in range(len(actions)):
            row.append([tdiag[l][q][j] - sum(ks[l][i][j]*generators[q][i][j] for i in range(rank)) for j in range(rank)])
        residuals.append(row)
    all_reached = set().union(*(set(col) for _, cols in raw_actions for col in cols))
    neighbors = sum(len(col) for _,cols in raw_actions for col in cols)
    return {'metric':g, 'generators':generators, 'projected':ks, 'residual_diagonal_blocks':residuals,
            'integer_actions':actions, 'action_denominators':dens,
            'amplitudes':amplitudes, 'rank':rank,
            'counts':{'selected_configurations':len(states), 'reached_configurations':len(all_reached),
                      'external_configurations':len(all_reached-set(states)), 'action_nonzeros':neighbors,
                      'fixed_N_dimension':comb(m,n), 'reduced_integer_entries':len(states)*rank,
                      'all_fixed_N_configurations_selected':len(states)==comb(m,n)}}


def comparison_envelope(moments):
    rank = moments['rank']
    g0 = moments['metric'][0][0]
    amplitudes = moments['amplitudes']
    aa = moments['generators']
    c = [[F(0) if i == j else ceil_grid(abs(aa[0][i][j]) +
          sum(a*abs(gen[i][j]) for a,gen in zip(amplitudes, aa[1:])))
          for j in range(rank)] for i in range(rank)]
    rr = moments['residual_diagonal_blocks']
    max_diag = [F(0)]*rank
    corners = list(cartesian_product(*[(-a,a) if a else (F(0),) for a in amplitudes]))
    for corner in corners:
        weights = (F(1),)+corner
        for j in range(rank):
            value = sum(weights[l]*weights[q]*rr[l][q][j] for l in range(len(weights)) for q in range(len(weights)))
            if value < 0:
                raise AssertionError('Negative exact residual square')
            max_diag[j] = max(max_diag[j], value)
    radii = [sqrt_up(x/g0) for x in max_diag]
    return c, radii, max_diag


def integrate_envelope(c, radii, horizon, order=24):
    """Bound integral r^T exp(C t)e0 dt by positive Taylor + scalar tail."""
    if horizon < 0 or type(order) is not int or order < 0:
        raise ValueError('Nonnegative horizon and order required')
    n = len(c)
    if len(radii) != n or any(len(row) != n or any(x < 0 for x in row) for row in c) or any(x < 0 for x in radii):
        raise ValueError('Nonnegative square comparison system required')
    if horizon == 0:
        return F(0), F(0)
    rate = max(sum(c[i][j] for i in range(n)) for j in range(n))
    x = rate*horizon
    ratio = x/F(order+3)
    if ratio >= 1:
        raise ValueError('Taylor order too small for geometric tail bound')
    term = [F(i == 0) for i in range(n)]
    total = F(0)
    for k in range(order+1):
        total += horizon*sum(a*b for a,b in zip(radii, term))/F(k+1)
        term = [horizon*sum(a*b for a,b in zip(row, term))/F(k+1) for row in c]
    # The omitted first scalar term is x^(order+1)/(order+2)!.
    scalar = F(1)
    for k in range(1, order+2):
        scalar *= x/F(k)
    tail = max(radii)*horizon*scalar/F(order+2)/(1-ratio)
    return total+tail, tail


def check(data, proposal, horizon=F(10), tolerance=F(1,200), order=24):
    started = time.monotonic()
    moments = exact_moments(data, proposal)
    c, radii, max_diag = comparison_envelope(moments)
    error, tail = integrate_envelope(c, radii, horizon, order)
    g = moments['metric']
    # Conventional uniform residual bound via Frobenius / Gershgorin metric.
    gmin = min(g[i][i]-sum(abs(v) for j,v in enumerate(g[i]) if i != j) for i in range(len(g)))
    naive = None if gmin <= 0 else horizon*sqrt_up(sum(max_diag)/gmin)
    energy = moments['projected'][0][0][0]/g[0][0]
    receipt = {'status':'accepted_uniform_control_reduction', 'fixture_sha256':digest(data),
        'proposal_sha256':digest(proposal), 'initial_state':'first exact rational subspace column, normalized',
        'ground_state_identity_claimed':False, 'controls':'arbitrary bounded measurable real functions within the declared amplitude box',
        'horizon_atomic_time':str(horizon), 'state_vector_error_bound':str(error), 'state_vector_error_float':float(error),
        'population_error_bound':str(min(F(2),2*error)), 'target_state_error':str(tolerance), 'target_met':error<=tolerance,
        'initial_energy_Ha':str(energy), 'initial_energy_float_Ha':float(energy),
        'envelope_tail':str(tail), 'comparison_matrix':[[str(x) for x in row] for row in c],
        'residual_column_radii':[str(x) for x in radii],
        'uniform_frobenius_error_float':None if naive is None else float(naive),
        'reduced_dimension':moments['rank'], **moments['counts'], 'exact_replay_seconds':time.monotonic()-started,
        'numeric_time_integration_included':False}
    return receipt, moments
