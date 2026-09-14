"""Reproducible finite-size diagnostic; exports no energy certificate."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import time

from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_general_schur import apply_columns, gram, fixture
from experiments.marginal_enlarged_schur import projected_action, reduced
from experiments.marginal_symbolic import add, adj, mono, scale, product, canonical
from experiments.marginal_symmetry_transfer import hopping_perturbation


def probe(modes, interaction):
    start = time.monotonic()
    m = modes // 2
    z = [{} for _ in range(m + 1)]
    for bits in range(1 << m):
        chosen = [i + m if bits & (1 << i) else i for i in range(m)]
        inversions = sum(chosen[i] > chosen[j] for i in range(m) for j in range(i + 1, m))
        z[bits.bit_count()][sum(1 << i for i in chosen)] = F((-1) ** inversions)
    assert gram(z, z) == [[F(comb(m, i)) if i == j else F(0) for j in range(m + 1)] for i in range(m + 1)]
    h0 = hopping_polynomial(modes, F(1, 5))
    assert not any(projected_action(h0, z)[2]), 'Reference embedding is not invariant'
    s = F(1, 100)
    h = add(h0, hopping_perturbation(modes, s)[0])
    for i in range(m):
        p = mono(((1, i), (0, (i + 1) % m + m)))
        h = add(h, scale(add(p, adj(p)), -s * F(m + i, m)))
    if interaction:
        p = mono(((1, 0), (1, 1), (0, m + 2), (0, m + 1)))
        n0, nr = mono(((1, 0), (0, 0))), mono(((1, m + 1), (0, m + 1)))
        h = add(h, scale(add(p, adj(p)), s / 2), scale(product(n0, nr), s / 7), scale(mono(((1, 2), (0, 2))), s / 11))
    h = canonical(h)
    if modes == 10:
        assert h == fixture(s, interaction), 'Ten-mode fixture mismatch'
    current = projected_action(h, z)[2]
    echelon = []
    ranks = []
    for degree in range(32):
        old_rank = len(echelon)
        for column in current:
            v = reduced(column, echelon)
            if v:
                if len(echelon) >= 128:
                    raise ValueError('Diagnostic rank budget exhausted')
                pivot = min(v)
                a = v[pivot]
                echelon.append((pivot, {s: c / a for s, c in v.items()}))
        ranks.append(len(echelon))
        if len(echelon) == old_rank:
            break
        current = apply_columns(h0, current)
    else:
        raise ValueError('Diagnostic power budget exhausted')
    r = [v for _, v in echelon]
    assert not any(reduced(v, echelon) for v in apply_columns(h0, r)), 'Closure residual'
    assert not any(a for row in gram(z, r) for a in row), 'Projection residual'
    if modes == 10:
        assert len(r) == (30 if interaction else 8)
    return {'modes': modes, 'particles': m, 'interaction': interaction, 'strength': str(s),
            'sector_dimension': comb(modes, m), 'reference_dimension': len(z),
            'coupling_ranks_by_power': ranks, 'retained_dimension': len(z) + len(r),
            'retained_support': len(set().union(*(set(v) for v in z + r))),
            'exact_closure_verified': True, 'exact_orthogonality_verified': True,
            'elapsed_seconds': time.monotonic() - start,
            'scope': 'Exact first-round dimension diagnostic only; no energy or scaling certificate'}


if __name__ == '__main__':
    for modes in (10, 12, 14, 16):
        for interaction in (False, True):
            out = Path(__file__).parent / ('dimension_m' + str(modes) + ('_mixed' if interaction else '_cycle') + '.json')
            if out.exists():
                raise ValueError('Preserve previous diagnostic')
            result = probe(modes, interaction)
            out.write_text(json.dumps(result, indent=2) + '\n')
            print(json.dumps(result), flush=True)
