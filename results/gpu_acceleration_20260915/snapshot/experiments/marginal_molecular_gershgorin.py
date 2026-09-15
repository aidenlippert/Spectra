"""Finite H4 certificates from automatic determinant selection and Gershgorin.

The explicit 70-state matrix is intentional. This tests reference discovery,
not implicit construction or asymptotic molecular scaling.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_symbolic import decode, hermitian
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_implicit_certificate import rational_text

ROOT = Path(__file__).resolve().parents[1]
SOURCES = {'rectangle': 'results/marginal_molecule/h4_compact_certificate.json',
           'square': 'results/marginal_molecule_stress/accepted_degree3_certificate.json'}


def matrix(certificate):
    if (type(certificate.get('modes')) is not int or certificate['modes'] != 8
            or type(certificate.get('particles')) is not int or certificate['particles'] != 4):
        raise ValueError('This finite diagnostic requires M8,N4')
    h = decode(certificate['hamiltonian'], 8, 4)
    if not hermitian(h) or any(sum(1 if c else -1 for c, _ in w) for w in h):
        raise ValueError('Real Hermitian number-conserving Hamiltonian required')
    states = [s for s in range(256) if s.bit_count() == 4]
    index = {s: i for i, s in enumerate(states)}
    a = [[F(0) for _ in states] for _ in states]
    for word, coefficient in h.items():
        for j, state in enumerate(states):
            result = apply_word(word, state)
            if result:
                destination, sign = result
                a[index[destination]][j] += coefficient * sign
    if any(a[i][j] != a[j][i] for i in range(70) for j in range(70)):
        raise ValueError('Nonsymmetric physical matrix')
    return states, a


def upper(a, witness):
    if type(witness) is not dict:
        raise ValueError('Missing integer upper witness')
    amplitudes = witness.get('amplitudes')
    if (type(amplitudes) is not list or len(amplitudes) != len(a)
            or any(type(x) is not int for x in amplitudes) or not any(amplitudes)):
        raise ValueError('Expected nonzero integer amplitudes in ascending bitstring order')
    norm = sum(x * x for x in amplitudes)
    energy = sum(x * y * a[i][j] for i, x in enumerate(amplitudes) for j, y in enumerate(amplitudes))
    return energy / norm


def row_bounds(a, q):
    return {i: a[i][i] - sum(abs(a[i][j]) for j in q if j != i) for i in q}


def prepare(a, retained):
    if (type(retained) is not list or not 1 <= len(retained) <= 32
            or any(type(i) is not int or not 0 <= i < len(a) for i in retained)
            or len(set(retained)) != len(retained)):
        raise ValueError('Expected one to 32 distinct retained indices')
    q = [i for i in range(len(a)) if i not in retained]
    gamma = min(row_bounds(a, q).values())
    block = [[a[i][j] for j in retained] for i in retained]
    leakage = [[sum(a[i][k] * a[k][j] for k in q) for j in retained] for i in retained]
    return block, leakage, gamma


def schur_pivots(data, lower):
    block, leakage, gamma = data
    lower = F(lower)
    if gamma <= lower:
        return None
    reduced = [[a - lower * (i == j) - leakage[i][j] / (gamma - lower)
                for j, a in enumerate(row)] for i, row in enumerate(block)]
    return ldl_pivots(reduced)


def suggest_lower(data, witness_upper):
    import numpy as np
    block, leakage, gamma = data
    a, k, c = np.array(block, dtype=float), np.array(leakage, dtype=float), float(gamma)
    def positive(b):
        return b < c and np.linalg.eigvalsh(a - b * np.eye(len(a)) - k / (c - b))[0] > 0
    hi = min(float(witness_upper), c)
    lo, step = hi - 1, 1
    for _ in range(32):
        if positive(lo): break
        step *= 2
        lo -= step
    else:
        raise ValueError('Numerical lower bracket exhausted')
    for _ in range(60):
        middle = (lo + hi) / 2
        if positive(middle): lo = middle
        else: hi = middle
    for exponent in range(10):
        lower = F(math.floor((lo - 10. ** (exponent - 10)) * 10**12), 10**12)
        if schur_pivots(data, lower) is not None and lower <= witness_upper:
            return lower
    raise ValueError('Exact Schur lower export failed')


def replay(certificate):
    if certificate.get('kind') != 'molecular_gershgorin_schur_v1':
        raise ValueError('Unsupported molecular Schur certificate')
    states, a = matrix(certificate)
    chosen = certificate.get('retained_states')
    if type(chosen) is not list or any(type(s) is not int or s not in states for s in chosen):
        raise ValueError('Invalid retained bitstrings')
    retained = [states.index(s) for s in chosen]
    data = prepare(a, retained)
    lower = F(certificate['lower'])
    checked = schur_pivots(data, lower)
    if checked is None:
        raise ValueError('Exact molecular Schur positivity failed')
    witness_upper = upper(a, certificate.get('independent_upper'))
    if witness_upper < lower:
        raise ValueError('Inconsistent molecular interval')
    return {'lower': rational_text(lower), 'upper': rational_text(witness_upper),
            'width': rational_text(witness_upper - lower), 'width_float': float(witness_upper - lower),
            'retained_dimension': len(retained), 'complement_dimension': len(a) - len(retained),
            'complement_lower': rational_text(data[2]), 'positive_denominator': rational_text(data[2] - lower),
            'schur_pivots': [rational_text(x) for x in checked],
            'scope': 'Exact rational 70-state Hamiltonian; explicit configuration matrix; no asymptotic scaling claim'}


def coupling_choice(a, p, q, data, lower):
    import numpy as np
    block, leakage, gamma = data
    effective = np.array(block, dtype=float) - np.array(leakage, dtype=float) / float(gamma - lower)
    _, vectors = np.linalg.eigh(effective)
    scores = [(sum(float(a[j][i]) * vectors[k, 0] for k, i in enumerate(p))) ** 2 for j in q]
    if not any(scores):
        rows = row_bounds(a, q)
        return min(q, key=lambda i: (rows[i], i))
    return q[max(range(len(q)), key=lambda k: (scores[k], -q[k]))]


def witness_choice(a, p, q, amplitudes):
    scores = {j: abs(sum(a[j][i] * amplitudes[i] for i in p)) for j in q}
    if any(scores.values()):
        return max(q, key=lambda i: (scores[i], -i))
    if any(amplitudes[j] for j in q):
        return max(q, key=lambda i: (abs(amplitudes[i]), -i))
    rows = row_bounds(a, q)
    return min(q, key=lambda i: (rows[i], i))


def run(name, target_width=F(1, 10**7), selector='row'):
    if name not in SOURCES:
        raise ValueError('Unknown H4 fixture')
    if selector not in ('row', 'coupling', 'witness'):
        raise ValueError('Unknown retained-state selector')
    out = ROOT / 'results/marginal_molecular_gershgorin' / (name + ('_' + selector if selector != 'row' else ''))
    if out.exists():
        raise ValueError('Preserve previous molecular run')
    out.mkdir(parents=True)
    started = time.monotonic()
    source = json.loads((ROOT / SOURCES[name]).read_text())
    states, a = matrix(source)
    witness_upper = upper(a, source['independent_upper'])
    p, q, history = [], list(range(70)), []
    best, lower = None, None
    while len(p) < 32:
        rows = row_bounds(a, q)
        if selector == 'coupling' and lower is not None:
            chosen = coupling_choice(a, p, q, data, lower)
        elif selector == 'witness' and lower is not None:
            chosen = witness_choice(a, p, q, source['independent_upper']['amplitudes'])
        else:
            chosen = min(q, key=lambda i: (rows[i], states[i]))
        p.append(chosen)
        q.remove(chosen)
        data = prepare(a, p)
        if data[2] <= witness_upper + F(1, 100):
            continue
        lower = suggest_lower(data, witness_upper)
        item = {'retained_dimension': len(p), 'gamma': str(data[2]),
                'lower': str(lower), 'width': str(witness_upper - lower), 'width_float': float(witness_upper - lower)}
        history.append(item)
        if best is None or lower > F(best['lower']):
            best = {'kind': 'molecular_gershgorin_schur_v1', 'modes': 8, 'particles': 4,
                    'hamiltonian': source['hamiltonian'], 'retained_states': [states[i] for i in p],
                    'lower': str(lower), 'independent_upper': {'amplitudes': source['independent_upper']['amplitudes']}}
        (out / 'history.json').write_text(json.dumps(history, indent=2) + '\n')
        print(json.dumps({k: item[k] for k in ('retained_dimension', 'width_float')}), flush=True)
        if witness_upper - lower <= target_width:
            break
    if best is None:
        raise ValueError('Retained-state budget exhausted before complement gap')
    receipt = replay(best)
    receipt.update(source=SOURCES[name], selector=selector, elapsed_seconds=time.monotonic() - started,
                   stopping_reason='target_width' if receipt['width_float'] <= float(target_width) else 'retained_dimension_budget',
                   target_width=str(target_width), final_greedy_dimension=len(p))
    (out / 'certificate.json').write_text(json.dumps(best, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('retained_dimension', 'width_float', 'elapsed_seconds', 'stopping_reason')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=tuple(SOURCES), default='rectangle')
    parser.add_argument('--verify')
    parser.add_argument('--selector', choices=('row', 'coupling', 'witness'), default='row')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    else:
        run(args.fixture, selector=args.selector)
