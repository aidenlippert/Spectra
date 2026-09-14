"""Exact support diagnostics, separate from numerical separator discovery."""
from fractions import Fraction as F
from itertools import combinations
import json
import sys
import time

from experiments.marginal_symbolic import canonical, mono
from research.certificate_scaling.commutator_dual_witness import moment_decode, psd
from research.joint_patterns_20260913.dual import integer_grams
from research.mechanism_transfer_20260913.ledger import modular_rank
from research.response_consistency_20260913.separator import ROOT, OUT


def run():
    start = time.monotonic()
    parent = json.loads((ROOT/'results/trace_pricing_20260913/full_dual/witness.json').read_text())
    y = moment_decode(parent['moments'], 12); records = []
    for support in combinations(range(6), 3):
        modes = [i for i in range(12) if i//2 in support]
        polys = [canonical(mono(tuple((0, i) for i in triple))) for triple in combinations(modes, 3)]
        matrices, _ = integer_grams(polys, [y], 'anticommutator')
        records.append({'spatial_support': list(support), **psd(matrices[0])})
    cert = json.loads((OUT/'separator.json').read_text()); pairs = list(combinations(range(12), 2))
    flattening = [[F(0) for _ in pairs] for _ in range(12)]; lookup = {pair: i for i, pair in enumerate(pairs)}
    for triple, z in zip(cert['triples'], cert['coefficients']):
        for position, mode in enumerate(triple):
            pair = tuple(i for i in triple if i != mode)
            flattening[mode][lookup[pair]] += (-1)**position*F(z, cert['denominator'])
    rank = modular_rank(flattening); occupied = len({i for triple in cert['triples'] for i in triple})
    if rank != occupied:
        raise ValueError('Exact support rank not established by the lower/upper sandwich')
    if any(n in sys.modules for n in ('numpy', 'scipy', 'cvxpy', 'pyscf')):
        raise AssertionError('Numerical imports during exact support audit')
    result = {'all_three_spatial_T1_blocks': records, 'selected_three_form_support_rank': rank,
        'support_rank_scope': 'This C cannot be expressed using fewer than eight one-particle modes under any invertible orbital change. Does not exclude another separator after an orbital rotation.',
        'many_body_states_enumerated': 0, 'numerical_packages_loaded': [], 'wall_seconds': time.monotonic()-start,
        'scope': 'No T1 violation exists on any three of the original six spatial orbitals, with arbitrary coefficients on that support.'}
    (OUT/'support_audit.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    run()
