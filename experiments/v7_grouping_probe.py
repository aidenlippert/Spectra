"""Preliminary exact norm grouping probe (not a full dynamics headroom claim)."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json
try:
    from experiments.pauli import commutator_i
    from experiments.certificates import PauliTerm, _anti, _sqrt_interval
except ModuleNotFoundError:
    from pauli import commutator_i
    from certificates import PauliTerm, _anti, _sqrt_interval
ROOT = Path(__file__).resolve().parents[1]

def chain_hamiltonian(n):
    h = {}
    for i in range(n - 1):
        for x, c in (('X', Fraction(1)), ('Z', Fraction(3, 2))):
            h[''.join(x if j in (i, i + 1) else 'I' for j in range(n))] = c
    for i in range(n): h[''.join('Z' if j == i else 'I' for j in range(n))] = Fraction(1, 2)
    return h

def residual_terms(n, order):
    """Exact coefficient of G^order O / order!, with no term truncation."""
    cur = {'X' + 'I' * (n - 1): Fraction(1)}
    for _ in range(order):
        nxt = defaultdict(Fraction)
        for op, coefficient in cur.items():
            for p, c in commutator_i(chain_hamiltonian(n), op).items(): nxt[p] += coefficient * c
        cur = nxt
    factorial = 1
    for k in range(2, order + 1): factorial *= k
    return [PauliTerm(p, c / factorial) for p, c in sorted(cur.items()) if p != 'I' * n and c]

def _pack(terms, order):
    groups = []; comparisons = 0
    for t in order:
        for g in groups:
            fits = True
            for u in g:
                comparisons += 1
                if not _anti(t.pauli, u.pauli): fits = False; break
            if fits: g.append(t); break
        else: groups.append([t])
    return groups, comparisons

def weighted_greedy(terms):
    remaining = list(terms); groups = []; comparisons = 0
    while remaining:
        seed = max(remaining, key=lambda t: abs(t.coeff)); remaining.remove(seed); group = [seed]
        while True:
            candidates = []
            for t in remaining:
                ok = True
                for u in group:
                    comparisons += 1
                    if not _anti(t.pauli, u.pauli): ok = False; break
                if ok: candidates.append(t)
            if not candidates: break
            t = max(candidates, key=lambda x: abs(x.coeff)); remaining.remove(t); group.append(t)
        groups.append(group)
    return groups, comparisons

def verify(groups, terms):
    expected = sorted((t.pauli, t.coeff) for t in terms); seen = []
    for g in groups:
        assert g and all(_anti(a.pauli, b.pauli) for a, b in combinations(g, 2))
        squared = sum((t.coeff * t.coeff for t in g), Fraction()); lo, hi = _sqrt_interval(squared, 12)
        assert lo * lo <= squared <= hi * hi; seen.extend((t.pauli, t.coeff) for t in g)
    assert sorted(seen) == expected

def measure(groups, terms, comparisons):
    verify(groups, terms)
    intervals = [_sqrt_interval(sum((t.coeff*t.coeff for t in g), Fraction()), 12) for g in groups]
    return {'groups': len(groups), 'norm_lower': str(sum((x[0] for x in intervals), Fraction())), 'norm_upper': str(sum((x[1] for x in intervals), Fraction())), 'norm_interval_width': str(sum((x[1]-x[0] for x in intervals), Fraction())), 'construction_comparisons': comparisons, 'checker_pair_checks': sum(len(g)*(len(g)-1)//2 for g in groups), 'sizes': sorted((len(g) for g in groups), reverse=True)}

def run():
    rows = []; failures = []
    for n in range(3, 9):
        for order in range(1, 6):
            terms = residual_terms(n, order)
            methods = {'firstfit': _pack(terms, terms), 'weight_greedy': weighted_greedy(terms), 'l1': ([[t] for t in terms], 0)}
            for name, (groups, comparisons) in methods.items():
                try: rows.append({'n': n, 'order': order, 'terms': len(terms), 'method': name, **measure(groups, terms, comparisons)})
                except AssertionError as exc: failures.append({'n': n, 'order': order, 'method': name, 'reason': str(exc)})
    return {'protocol': {'instances': 30, 'source': 'exact G^m O/m! from local XX+ZZ+Z chain', 'coefficient_type': 'Fraction', 'claim_scope': 'norm grouping only'}, 'rows': rows, 'failures': failures}

if __name__ == '__main__':
    out = run(); path = ROOT/'results/v7/grouping_probe.json'; path.parent.mkdir(exist_ok=True); path.write_text(json.dumps(out, indent=2)); print(json.dumps(out, indent=2))
