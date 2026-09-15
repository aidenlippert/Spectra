"""Exact, non-enumerating upper bounds for disjoint orbital rotations.

Only one-particle projectors and CAR Wick contractions are constructed. The
unnormalized occupied columns fix the norm, so a compact recipe is not confused
with an implicitly enumerated state vector. This is a supporting control, not
the full response-dressed upper construction.
"""
from fractions import Fraction as F
import hashlib
import json
import time

from experiments.marginal_symbolic import decode, hermitian


def digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def projector(spatial, occupied, tangents):
    if type(spatial) is not int or spatial < 2 or spatial % 2:
        raise ValueError('An even number of spatial orbitals is required')
    if type(occupied) is not int or not 0 <= occupied <= spatial:
        raise ValueError('Invalid occupied orbital count')
    if len(tangents) != spatial // 2:
        raise ValueError('One tangent for each disjoint spatial pair is required')
    ts = [F(t) for t in tangents]
    if any(max(abs(t.numerator).bit_length(), t.denominator.bit_length()) > 256 for t in ts):
        raise ValueError('Tangent bit budget exceeded')
    r = [[F(0) for _ in range(spatial)] for _ in range(spatial)]
    norm = F(1)
    for i, t in enumerate(ts):
        j = i + spatial // 2
        for col, vector in ((i, (F(1), t)), (j, (-t, F(1)))):
            if col >= occupied:
                continue
            norm *= 1 + t*t
            for a, va in zip((i, j), vector):
                for b, vb in zip((i, j), vector):
                    r[a][b] += va*vb/(1+t*t)
    if sum(r[i][i] for i in range(spatial)) != occupied:
        raise AssertionError('Projector particle trace failed')
    if any(sum(r[i][k]*r[k][j] for k in range(spatial)) != r[i][j]
           for i in range(spatial) for j in range(spatial)):
        raise AssertionError('Occupied projector is not idempotent')
    return r, norm


def wick(word, r):
    """Expectation with r[p][q] = <a_q^dagger a_p>, including nonnormal words."""
    if not word:
        return F(1)
    if len(word) % 2:
        return F(0)
    ca, ia = word[0]
    value = F(0)
    for j in range(1, len(word)):
        cb, ib = word[j]
        if ca == cb:
            continue
        contraction = r[ib][ia] if ca else F(ia == ib)-r[ia][ib]
        if contraction:
            value += (-1)**(j-1)*contraction*wick(word[1:j]+word[j+1:], r)
    return value


def evaluate(data, alpha_t, beta_t, na, nb):
    m = data['modes']
    if type(m) is not int or m < 4 or m % 4 or type(data['particles']) is not int:
        raise ValueError('Invalid disjoint-pair model')
    if type(na) is not int or type(nb) is not int or na+nb != data['particles']:
        raise ValueError('Particle sector mismatch')
    s = m//2
    ra, norma = projector(s, na, alpha_t)
    rb, normb = projector(s, nb, beta_t)
    r = [[F(0) for _ in range(m)] for _ in range(m)]
    for spin, rs in enumerate((ra, rb)):
        for i in range(s):
            for j in range(s):
                r[2*i+spin][2*j+spin] = rs[i][j]
    h = decode(data['hamiltonian'], m, 4)
    if not hermitian(h) or any(sum(2*c-1 for c, _ in w) for w in h):
        raise ValueError('A Hermitian, particle-conserving Hamiltonian is required')
    e = sum((c*wick(w, r) for w, c in h.items()), F(0))
    norm = norma*normb
    s2 = F((na-nb)**2, 4) + F(na+nb, 2) - sum(ra[i][j]*rb[j][i] for i in range(s) for j in range(s))
    return {'energy': e, 'norm': norm, 'numerator': e*norm, 'S2': s2,
            'hamiltonian_terms': len(h), 'projector_entries': m*m,
            'many_body_states_enumerated': 0}


def make_certificate(data, alpha_t, beta_t, na, nb, target_spin=None):
    a, b = list(map(F, alpha_t)), list(map(F, beta_t))
    result = evaluate(data, a, b, na, nb)
    cert = {'kind': 'disjoint_slater_upper_v1', 'fixture_sha256': digest(data),
            'n_alpha': na, 'n_beta': nb, 'alpha_tangents': list(map(str, a)),
            'beta_tangents': list(map(str, b)), 'target_spin': target_spin,
            'upper_Ha': str(result['energy']), 'norm_squared': str(result['norm'])}
    check(data, cert)
    return cert


def check(data, cert):
    start = time.monotonic()
    keys = {'kind', 'fixture_sha256', 'n_alpha', 'n_beta', 'alpha_tangents',
            'beta_tangents', 'target_spin', 'upper_Ha', 'norm_squared'}
    if set(cert) != keys or cert['kind'] != 'disjoint_slater_upper_v1':
        raise ValueError('Unknown compact upper schema')
    if cert['fixture_sha256'] != digest(data):
        raise ValueError('Compact upper input binding failed')
    if any(type(x) is not str for x in [cert['upper_Ha'], cert['norm_squared'],
                                      *cert['alpha_tangents'], *cert['beta_tangents']]):
        raise ValueError('Exact rational strings are required')
    result = evaluate(data, cert['alpha_tangents'], cert['beta_tangents'],
                      cert['n_alpha'], cert['n_beta'])
    if F(cert['upper_Ha']) != result['energy'] or F(cert['norm_squared']) != result['norm'] or result['norm'] <= 0:
        raise ValueError('Upper energy or norm did not reproduce')
    spin = cert['target_spin']
    if spin is not None:
        if type(spin) is not int or spin not in (0, 1) or result['S2'] != spin*(spin+1):
            raise ValueError('Requested pure-spin condition failed')
        # The supported sufficient spin certificate uses a common orbital
        # frame and nested occupied alpha/beta spaces. <S2> alone can be
        # insufficient for an intermediate spin in a mixture.
        if cert['n_alpha']-cert['n_beta'] != 2*spin or list(map(F, cert['alpha_tangents'])) != list(map(F, cert['beta_tangents'])):
            raise ValueError('Pure-spin replay requires nested common orbitals')
    return {'upper_Ha': str(result['energy']), 'norm_squared': str(result['norm']),
            'numerator_Ha': str(result['numerator']), 'S2': str(result['S2']),
            'spin_certificate': 'common orbitals, nested occupations' if spin is not None else None,
            'hamiltonian_terms': result['hamiltonian_terms'], 'projector_entries': result['projector_entries'],
            'many_body_states_enumerated': 0, 'replay_seconds': time.monotonic()-start}
