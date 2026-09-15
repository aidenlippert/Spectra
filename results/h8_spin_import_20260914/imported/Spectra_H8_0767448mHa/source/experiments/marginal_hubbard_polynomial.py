"""Direct local-atom Hubbard-chain certificates for the existing exact checker.

No configuration list, numerical optimizer, or full-population ideal lift is
used. Expanded positive atoms number O(sites**2); their degree is at most four.
"""
from fractions import Fraction as F
from math import lcm
import json
from pathlib import Path

from experiments.marginal_symbolic import mono, product, scale, add, encode
from experiments.marginal_polynomial_metric import replay


def build(sites, U, t, slack=F(1, 1000)):
    """Certify QHQ >= U-4*t*(sites-1)-slack for the neutral spin-zero chain."""
    if type(sites) is not int or not 2 <= sites <= 32 or sites % 2:
        raise ValueError('An even number of sites in 2..32 is required')
    if any(type(x) not in (int, str, F) for x in (U, t, slack)):
        raise ValueError('Exact rational parameters required')
    U, t, slack = F(U), F(t), F(slack)
    if U < 0 or t < 0 or slack <= 0:
        raise ValueError('Nonnegative U,t and positive slack required')
    edges = sites-1
    gamma = U-4*t*edges-slack
    denominator = lcm(2*U.denominator, 2*t.denominator,
                      2*gamma.denominator, slack.denominator)
    if denominator > 10**16:
        raise ValueError('Positivity denominator exceeds existing verifier budget')
    n = [mono(((1, i), (0, i))) for i in range(2*sites)]
    h = add(*(scale(product(n[2*i], n[2*i+1]), U) for i in range(sites)),
            *(mono(((1, 2*i+s), (0, 2*(i+1)+s)), -t)
              for i in range(edges) for s in (0, 1)),
            *(mono(((1, 2*(i+1)+s), (0, 2*i+s)), -t)
              for i in range(edges) for s in (0, 1)))
    positive = {}

    def indicator(occupied, empty, value):
        if occupied & empty or not value:
            return
        key = (occupied | empty, occupied)
        positive[key] = positive.get(key, F(0))+value

    # K_D=U D(D-1)+(2*t*edges+slack)(D-1)+slack+t*sum(S).
    # S=(1-If-Ib)(D+1)+If(1-ni_bar+nj_bar)+Ib(1-nj_bar+ni_bar).
    # Equal spin occupations give 1-If-Ib=ni*nj+(1-ni)*(1-nj).
    for i in range(edges):
        for spin in (0, 1):
            ni, nj = 1 << (2*i+spin), 1 << (2*(i+1)+spin)
            bi, bj = 1 << (2*i+1-spin), 1 << (2*(i+1)+1-spin)
            for occupied, empty in [(ni | nj, 0), (0, ni | nj)]:
                indicator(occupied, empty, t)
                for k in range(sites):
                    indicator(occupied | (3 << (2*k)), empty, t)
            indicator(nj, ni | bi, t)
            indicator(nj | bj, ni, t)
            indicator(ni, nj | bj, t)
            indicator(ni | bi, nj, t)

    # v=D-T/2 with T=(Nalpha-target)+(Nbeta-target).
    # K_v=K_D-T/2*(U*D-gamma-t*sum_directed I).
    ideal = {0: gamma/2}
    for i in range(sites):
        ideal[3 << (2*i)] = -U/2
    for i in range(edges):
        for spin in (0, 1):
            ni, nj = 1 << (2*i+spin), 1 << (2*(i+1)+spin)
            for mask, value in [(ni, t/2), (nj, t/2), (ni | nj, -t)]:
                ideal[mask] = ideal.get(mask, F(0))+value
    multipliers = [{'mask': mask, 'coefficient': int(value*denominator)}
                   for mask, value in sorted(ideal.items()) if value]
    localizers = [{'required': 0, 'occupied': 0,
                   'weight': int((2*t*edges+slack)*denominator)}]
    if U:
        localizers += [{'required': 3 << (2*i), 'occupied': 3 << (2*i),
                        'weight': int(U*denominator)} for i in range(sites)]
    return {
        'kind': 'joint_polynomial_metric_gap_v1', 'modes': 2*sites,
        'particles': sites, 'hamiltonian': encode(h), 'target_lower': str(gamma),
        'polynomial_metric': {'denominator': 2, 'terms': [
            {'powers': [2 if i == j else 0 for i in range(sites)], 'coefficient': 1}
            for j in range(sites)]},
        'weight_proof': {'denominator': 2, 'bound': 2,
            'positive_indicators': [],
            'charge_indicators': [{'required': 0, 'occupied': 0, 'weight': 2}],
            'number_multipliers': [[{'mask': 0, 'coefficient': -1}] for _ in (0, 1)]},
        'numerator_proof': {'denominator': denominator,
            'bound': int(slack*denominator),
            'positive_indicators': [{'required': required, 'occupied': occupied,
                                     'weight': int(value*denominator)}
                                    for (required, occupied), value in sorted(positive.items())],
            'charge_indicators': localizers,
            'number_multipliers': [list(multipliers), list(multipliers)]}}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('sites', type=int)
    parser.add_argument('--U', default='4')
    parser.add_argument('--t', default='1/3')
    parser.add_argument('--slack', default='1/1000')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    certificate = build(args.sites, args.U, args.t, args.slack)
    receipt = replay(certificate)
    for part in ('weight_positivity', 'numerator_positivity'):
        if F(receipt[part]['residual_l1']):
            raise ValueError('Analytic construction must have zero exact residual')
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out/'certificate.json').write_text(json.dumps(certificate, indent=2)+'\n')
    (args.out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))
