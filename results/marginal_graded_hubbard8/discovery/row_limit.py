"""Exact finite-chain D-metric row diagnostic; not a spectrum solver."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_polynomial_metric import JointPolynomial


def penalty(left, right, D):
    total = 0
    before = int(left == 3)+int(right == 3)
    for bit in (1, 2):
        if bool(left&bit) != bool(right&bit):
            after = int((left^bit) == 3)+int((right^bit) == 3)
            total += D+after-before
    return total


def limits(sites, U=4, t=1):
    population = sites//2; rows = []
    for D in range(1, population+1):
        current = {(0, 0, 0, -1): (0, ())}; visited = 0
        for position in range(sites):
            following = {}; remaining = sites-position-1
            for (a, b, d, last), (score, word) in current.items():
                visited += 1
                for symbol in range(4):
                    aa, bb, dd = a+bool(symbol&1), b+bool(symbol&2), d+(symbol == 3)
                    if not aa <= population <= aa+remaining or not bb <= population <= bb+remaining or not dd <= D <= dd+remaining:
                        continue
                    key = (aa, bb, dd, symbol)
                    value = score+(penalty(last, symbol, D) if last >= 0 else 0)
                    if key not in following or value > following[key][0]:
                        following[key] = (value, word+(symbol,))
            current = following
        score, word = max(current.values())
        lower = F(U)*D-F(t)*score/D
        rows.append({'D': D, 'maximum_weighted_hopping_sum': score, 'row_lower': str(lower),
                     'one_maximizer': list(word), 'DP_states_processed': visited})
    return rows


if __name__ == '__main__':
    root = Path('results/marginal_graded_hubbard8')
    rows = limits(8)
    worst = min(rows, key=lambda r: F(r['row_lower']))
    ring = JointPolynomial(build(8, 4, 1)); weight, numerator, cost = ring.compile(F(0))
    state = sum(symbol << (2*i) for i, symbol in enumerate(worst['one_maximizer']))
    evaluate = lambda p, s: sum(v for m, v in p.items() if s&m == m)
    w = F(evaluate(weight, state), cost['metric_scale'])
    energy = F(evaluate(numerator, state), cost['numerator_scale'])/w
    assert energy == F(worst['row_lower'])
    attempted = energy+F(1, 100)
    _, attempted_numerator, attempted_cost = ring.compile(attempted)
    violated = F(evaluate(attempted_numerator, state), attempted_cost['numerator_scale'])
    assert violated < 0
    small = JointPolynomial(build(4, 4, 1)); sw, sk, sc = small.compile(F(0)); brute = {}
    for a in combinations(range(4), 2):
        for b in combinations(range(4), 2):
            s = sum(1 << (2*i) for i in a)+sum(1 << (2*i+1) for i in b)
            D = len(set(a)&set(b))
            if not D:
                continue
            value = F(evaluate(sk, s), sc['numerator_scale'])/F(evaluate(sw, s), sc['metric_scale'])
            brute[D] = min(brute.get(D, value), value)
    assert brute == {r['D']: F(r['row_lower']) for r in limits(4)}
    receipt = {'sites': 8, 'U': '4', 't': '1', 'metric': 'D', 'sector': 'Nalpha=Nbeta=4, D>=1',
               'rows_by_D': rows, 'exact_fixed_metric_weighted_row_limit': str(energy),
               'critical_state': state, 'attempted_gamma': str(attempted),
               'exact_negative_numerator_at_critical_state': str(violated),
               'crosscheck': 'Four-site fixed-spin polynomial rows exhaustively match DP; eight-site critical row matches exact compiled polynomial.',
               'scope': 'Exact minimum of D-weighted Gershgorin rows for this Hubbard chain, via an integer local-state DP. Not the actual complement eigenvalue. A negative physical numerator rules out all nonnegative numerator decompositions at the attempted gamma for this fixed metric (and its positive scalar multiples). Other metrics and operator certificates remain possible.'}
    (root/'fixed_metric_row_limit.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt))
