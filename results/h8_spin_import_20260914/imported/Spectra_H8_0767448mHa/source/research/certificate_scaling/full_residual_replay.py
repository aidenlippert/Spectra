"""Exact full-certificate residual regrouping; standard-library replay only.

This research checker supplements, rather than changes, the production checker.
For R=cI+dGamma(A)+T and alpha=tr(A)/M, the fixed-N bound is
R >= c+N*alpha-min(N,M-N)*||A-alpha I|| - ||T||.
Use rational row sums for the matrix norm and CAR coefficient l1 for T.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import (
    add, decode, expand_squares, expand_orbits, hermitian, mono,
    number_shift, product, scale, verify,
)


def regroup(residual, modes, particles):
    if not hermitian(residual):
        raise ValueError('Residual must be Hermitian')
    one = [[F(0) for _ in range(modes)] for _ in range(modes)]
    tail = F(0)
    for w, c in residual.items():
        if not w:
            continue
        if len(w) == 2 and w[0][0] == 1 and w[1][0] == 0:
            one[w[0][1]][w[1][1]] += c
        else:
            tail += abs(c)
    alpha = sum(one[i][i] for i in range(modes)) / modes
    row_norm = max(sum(abs(one[i][j] - (alpha if i == j else 0))
                       for j in range(modes)) for i in range(modes))
    constant = residual.get((), F(0)) + particles * alpha
    allowance = min(particles, modes-particles) * row_norm + tail
    return {'constant_on_sector': constant, 'centered_onebody_row_norm': row_norm,
            'higher_body_l1': tail, 'residual_lower': constant - allowance}


def replay(path):
    start = time.monotonic()
    raw = path.read_bytes()
    cert = json.loads(raw)
    baseline = verify(cert)
    modes, particles = cert['modes'], cert['particles']
    degree = cert.get('operator_degree', 3)
    h = decode(cert['hamiltonian'], modes, 4)
    x = (decode(cert['number_multiplier'], modes, 2*degree-2)
         if 'number_multiplier' in cert else expand_orbits(cert['multiplier_orbits'], modes))
    squares, _ = expand_squares(cert['blocks'], cert['denominator'], modes, degree)
    r = add(h, mono((), -F(cert['b'])), scale(squares, -1),
            scale(product(number_shift(modes, particles), x), -1))
    grouped = regroup(r, modes, particles)
    lower = max(F(baseline['lower']), F(cert['b']) + grouped['residual_lower'])
    return {'certificate': str(path), 'sha256': hashlib.sha256(raw).hexdigest(),
            'baseline_lower': baseline['lower'], 'structured_lower': str(lower),
            'improvement': str(lower-F(baseline['lower'])),
            'improvement_float': float(lower-F(baseline['lower'])),
            'details': {k:str(v) for k,v in grouped.items()},
            'seconds': time.monotonic()-start, 'full_sector_enumeration': False,
            'status': 'exact_rational_research_bound; production_checker_unchanged'}


def self_test():
    n0 = ((1,0),(0,0)); n1 = ((1,1),(0,1))
    assert regroup({():F(3),n0:F(2),n1:F(2)},2,1)['residual_lower'] == 5
    hopping = {((1,0),(0,1)):F(2),((1,1),(0,0)):F(2)}
    assert regroup(hopping,2,1)['residual_lower'] == -2
    assert regroup(hopping,2,0)['residual_lower'] == 0
    assert regroup(hopping,2,2)['residual_lower'] == 0
    assert regroup({n0:F(1),n1:F(-1)},2,1)['residual_lower'] == -1


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('certificates', nargs='*', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    self_test()
    result = {'self_tests':'passed', 'receipts':[replay(q) for q in a.certificates]}
    text = json.dumps(result, indent=2)+'\n'
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(text)
    print(text)
