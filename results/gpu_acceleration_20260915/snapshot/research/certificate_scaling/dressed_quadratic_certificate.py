"""Proposal and independent rational replay for a restricted hidden-free family.

Only proposal imports NumPy. Exact replay regenerates H and verifies both a
standard SOS lower certificate and a rational dressed Slater upper witness.
"""
from fractions import Fraction as F
from experiments.marginal_symbolic import encode, decode, canonical, verify
from research.certificate_scaling.dressed_quadratic_recognition import dress_polynomial, recognize


def rational(x):
    if type(x) not in (int, str, F):
        raise ValueError('Expected an exact rational')
    return F(x)


def quadratic_matrix(q, modes):
    if type(modes) is not int or modes < 1 or not isinstance(q, dict):
        raise ValueError('Invalid quadratic input')
    A = [[F(0) for _ in range(modes)] for _ in range(modes)]
    c0 = F(0)
    for w, value in q.items():
        c = rational(value)
        if not w:
            c0 += c
            continue
        if (type(w) is not tuple or len(w) != 2 or
            any(type(t) is not tuple or len(t) != 2 or type(t[0]) is not int or
                type(t[1]) is not int or not 0 <= t[1] < modes for t in w) or
            w[0][0] != 1 or w[1][0] != 0):
            raise ValueError('Expected canonical number-conserving quadratic terms')
        A[w[0][1]][w[1][1]] += c
    if any(A[i][j] != A[j][i] for i in range(modes) for j in range(modes)):
        raise ValueError('Non-Hermitian quadratic input')
    return c0, A


def rational_slater_upper(q, C, modes, particles):
    c0, A = quadratic_matrix(q, modes)
    if type(particles) is not int or not 0 <= particles <= modes:
        raise ValueError('Invalid particle sector')
    if not isinstance(C, list) or len(C) != modes or any(not isinstance(r, list) or len(r) != particles for r in C):
        raise ValueError('Occupied columns have wrong dimensions')
    C = [[rational(x) for x in row] for row in C]
    # Exact Gauss-Jordan inverse of the occupied Gram; no sector enumeration.
    n = particles
    gram = [[sum(C[k][i]*C[k][j] for k in range(modes)) for j in range(n)] for i in range(n)]
    aug = [row[:] + [F(i == j) for j in range(n)] for i, row in enumerate(gram)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if aug[i][j]), None)
        if pivot is None:
            raise ValueError('Occupied columns are rank deficient')
        aug[j], aug[pivot] = aug[pivot], aug[j]
        d = aug[j][j]
        aug[j] = [v/d for v in aug[j]]
        for i in range(n):
            if i != j:
                d = aug[i][j]
                aug[i] = [v-d*w for v,w in zip(aug[i], aug[j])]
    inv = [row[n:] for row in aug]
    # Trace(C^T A C G^-1) using O(M^2 N + M N^2 + N^3) operations.
    AC = [[sum(A[i][j]*C[j][k] for j in range(modes)) for k in range(n)] for i in range(modes)]
    B = [[sum(C[k][i]*AC[k][j] for k in range(modes)) for j in range(n)] for i in range(n)]
    return c0 + sum(B[i][j]*inv[j][i] for i in range(n) for j in range(n))


def propose(h, modes, particles, denominator=10**8):
    import numpy as np
    if type(particles) is not int or not 0 <= particles <= modes:
        raise ValueError('Invalid particle sector')
    if type(denominator) is not int or denominator < 1:
        raise ValueError('Invalid denominator')
    recognition = recognize(h, modes)
    if not recognition.get('accepted'):
        raise ValueError('Structural recognition refused: '+recognition.get('reason', 'unknown'))
    q = recognition['recovered_quadratic']
    matching = recognition['matching']
    c0, A = quadratic_matrix(q, modes)
    vals, vecs = np.linalg.eigh(np.array(A, dtype=float))
    if not np.isfinite(vals).all() or not np.isfinite(vecs).all():
        raise ValueError('Non-finite eigenproposal')
    target_mu = vals[0]-1 if particles == 0 else vals[-1]+1 if particles == modes else (vals[particles-1]+vals[particles])/2
    mu = F(round(float(target_mu)*denominator), denominator)
    blocks = []
    qtrace = F(0)
    for k, val in enumerate(vals):
        delta = float(val)-float(mu)
        row = [round(float(x)*abs(delta)**0.5*denominator) for x in vecs[:,k]]
        if not any(row):
            continue
        # p=a yields p†p=a†a; hole p=a† yields aa†=1-a†a.
        creation = int(delta < 0)
        if creation:
            qtrace += F(sum(x*x for x in row), denominator**2)
        p = {((creation, i),): F(x) for i,x in enumerate(row) if x}
        dressed = canonical(dress_polynomial(p, matching, modes))
        words = sorted(dressed, key=lambda w:(len(w),w))
        if any(dressed[w].denominator != 1 for w in words):
            raise AssertionError('Integer CZ factor lost integrality')
        blocks.append({'words':[[list(t) for t in w] for w in words],
                       'factor':[[int(dressed[w]) for w in words]]})
    certificate = {'modes':modes, 'particles':particles, 'hamiltonian':encode(h),
                   'b':str(c0+mu*particles-qtrace), 'number_multiplier':encode({():mu}),
                   'denominator':denominator, 'blocks':blocks}
    C = [[str(F(round(float(vecs[i,j])*denominator), denominator)) for j in range(particles)] for i in range(modes)]
    payload = {'certificate':certificate, 'quadratic':encode(q),
               'matching':[list(p) for p in matching], 'occupied_columns':C,
               'recognition_work':{k:v for k,v in recognition.items() if k.endswith('work')},
               'scope':'Restricted matching-CZ-dressed quadratic class; no general molecular accuracy theorem.'}
    replay(payload)
    return payload


def replay(payload):
    cert = payload['certificate']
    m, n = cert['modes'], cert['particles']
    # The production checker owns lower-bound validation, including sector/type gates.
    lower_receipt = verify(cert)
    h = decode(cert['hamiltonian'], m, 4)
    q = decode(payload['quadratic'], m, 2)
    quadratic_matrix(q, m)
    raw_matching = payload['matching']
    if not isinstance(raw_matching, list) or any(not isinstance(p, list) for p in raw_matching):
        raise ValueError('Invalid matching encoding')
    matching = tuple(tuple(p) for p in raw_matching)
    if canonical(dress_polynomial(q, matching, m)) != h:
        raise ValueError('Dressed quadratic does not regenerate original H')
    hi = rational_slater_upper(q, payload['occupied_columns'], m, n)
    lo = F(lower_receipt['lower'])
    if hi < lo:
        raise AssertionError('Inconsistent exact interval')
    return {'lower':str(lo), 'upper':str(hi), 'width':str(hi-lo),
            'lower_float':float(lo), 'upper_float':float(hi), 'width_float':float(hi-lo),
            'passes_0_0016_Ha':hi-lo <= F(16,10000), 'lower_replay':lower_receipt,
            'upper_method':'exact rational dressed Slater projector', 'modes':m, 'particles':n,
            'scope':'Both bounds on original rational H; hidden-free structural control only.'}

if __name__ == '__main__':
    import argparse, json, time
    from pathlib import Path
    parser=argparse.ArgumentParser()
    parser.add_argument('--witness',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();start=time.monotonic()
    result=replay(json.loads(args.witness.read_text()))
    result['replay_seconds']=time.monotonic()-start
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
