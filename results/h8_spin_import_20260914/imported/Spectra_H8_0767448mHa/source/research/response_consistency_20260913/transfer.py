"""One independent molecular counterexample test; no T1 energy enrichment."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

from experiments.marginal_symbolic import add, canonical, mono, scale
from research.certificate_scaling.commutator_dual_witness import moment_decode, evaluate, seed
from research.joint_patterns_20260913.core import anticommutator
from research.joint_patterns_20260913.dual import affine_round, integer_grams
from research.molecular_collective_20260913.core import digest, extract, propose_tail
from research.spin_enrichment_20260913.full_dual import full_groups, check as full_check
from research.response_consistency_20260913.separator import ROOT, OUT, operator, domain

DEST = OUT/'fresh_h6_1p6'


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def discover():
    start = time.monotonic()
    import numpy as np
    from scipy.linalg import eigh
    from research.mechanism_transfer_20260913.discovery import Model
    data = json.loads((DEST/'fixture.json').read_text()); p = extract(data)
    tail, tail_discovery = propose_tail(data, p, 10); save(DEST/'tail.json', tail)
    model = Model(data, tail); span = model.full_span()
    if len(span) != 492:
        raise ValueError('Fresh case does not have the prescribed completed frame')
    primal, raw = model.solve(span, DEST/'old_full_solve', 300-(time.monotonic()-start), solver_cap=45.)
    rounded, rank = affine_round(raw, 12, 6); words = domain(12)
    trace = {w: seed(w, 12, 6) for w in words}
    base, added = full_groups(model.p, tail, model.symmetry['parity_masks']); thresholds = []
    for group in base+added:
        matrices, den = integer_grams(group['polynomials'], [rounded, trace], group['kind'])
        candidate, anchor = [np.array([[float(F(x, den)) for x in row] for row in a]) for a in matrices]
        value = float(eigh(candidate, anchor, eigvals_only=True, subset_by_index=[0, 0])[0])
        thresholds.append(max(0., -value/(1-value)) if value < 0 else 0.)
    mixture = F(max(10**6, math.ceil(1.1*max(thresholds)*10**12)), 10**12)
    save(DEST/'mixture_proposal.json', {'trace_mixture': str(mixture), 'maximum_required_fraction_float': max(thresholds),
        'thresholds': thresholds, 'affine_rank': rank, 'seconds': time.monotonic()-start})
    if mixture > F(1, 100):
        raise ValueError('Fresh full-family repair exceeds the declared one-percent cap')
    y = {w: (1-mixture)*rounded.get(w, F(0))+mixture*trace[w] for w in words}
    witness = {'kind': 'spin_completion_full_dual_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'parity_masks': model.symmetry['parity_masks'], 'trace_mixture': str(mixture),
        'moments': [{'word': w, 'value': str(v)} for w, v in sorted(y.items(), key=lambda x: (len(x[0]), x[0])) if v or not w]}
    save(DEST/'witness.json', witness)
    frozen = json.loads((OUT/'separator.json').read_text()); _, P = operator(frozen)
    fixed_value = evaluate(P, y); selected = None; support_records = []
    for count in (3, 4):
        candidates = []
        for support in combinations(range(6), count):
            triples = list(combinations([i for i in range(12) if i//2 in support], 3))
            polys = [canonical(mono(tuple((0, i) for i in t))) for t in triples]
            matrices, den = integer_grams(polys, [y], 'anticommutator')
            G = np.array([[float(F(x, den)) for x in row] for row in matrices[0]])
            vals, vecs = np.linalg.eigh(G); candidates.append((float(vals[0]), support, triples, polys, vecs[:, 0]))
        value, support, triples, polys, vector = min(candidates, key=lambda x: x[0])
        record = {'spatial_count': count, 'subsets_tested': len(candidates), 'minimum_eigenvalue': value,
                  'best_support': support, 'rounding_attempts': []}; support_records.append(record)
        if value >= -1e-9:
            continue
        for den in (10, 100, 1000, 10**6):
            coefficients = np.rint(vector*den).astype(int)
            C = add(*(scale(poly, F(int(z), den)) for poly, z in zip(polys, coefficients) if z))
            Q = anticommutator(C, C); exact = evaluate(Q, y)
            attempt = {'denominator': den, 'exact_moment': str(exact),
                'triples': [list(t) for t, z in zip(triples, coefficients) if z],
                'coefficients': [int(z) for z in coefficients if z], 'accepted': exact < 0}
            record['rounding_attempts'].append(attempt)
            if exact < 0:
                selected = {**frozen, 'parent_witness_sha256': hashlib.sha256((DEST/'witness.json').read_bytes()).hexdigest(),
                    'support_spatial': list(support), **{k: attempt[k] for k in ('denominator', 'exact_moment', 'triples', 'coefficients')},
                    'moment_float': float(exact), 'quartic_terms': len(Q), 'maximum_degree': max(map(len, Q))}
                save(DEST/'separator.json', selected); break
        if selected:
            break
    save(DEST/'discovery.json', {'geometry_spacing_Angstrom': '1.6', 'FCI_used': False,
        'model_construction': model.construction, 'old_full_solve': {k: primal[k] for k in ('wall_seconds', 'solve_seconds', 'status', 'directions')},
        'tail_discovery': tail_discovery, 'trace_mixture': str(mixture),
        'frozen_C_exact_moment': str(fixed_value), 'frozen_C_rejects': fixed_value < 0,
        'parametric_T1_search': support_records, 'parametric_separator_found': selected is not None,
        'wall_seconds': time.monotonic()-start, 'many_body_states_enumerated': 0,
        'scope': 'One new fixed-full-family solve proposes a witness. Exact acceptance is a separate process; no T1 energy optimization.'})


def replay():
    start = time.monotonic()
    from research.response_consistency_20260913.separator import check as separator_check
    data = json.loads((DEST/'fixture.json').read_text()); tail = json.loads((DEST/'tail.json').read_text())
    blob = (DEST/'witness.json').read_bytes(); parent = json.loads(blob)
    accepted = full_check(data, tail, parent)
    frozen = json.loads((OUT/'separator.json').read_text()); _, P = operator(frozen)
    y = moment_decode(parent['moments'], 12); fixed_value = evaluate(P, y)
    selected_path = DEST/'separator.json'
    separator = separator_check(json.loads(selected_path.read_text()), blob) if selected_path.exists() else None
    if separator:
        separator['scope'] = 'Exact separation of the fresh hash-bound functional; full old-family feasibility is checked in this same receipt.'
    if any(n in sys.modules for n in ('numpy', 'scipy', 'cvxpy', 'pyscf')):
        raise AssertionError('Numerical imports in the fresh molecular replay')
    save(DEST/'exact_receipt.json', {'full_old_family': accepted,
        'fixture_sha256': hashlib.sha256((DEST/'fixture.json').read_bytes()).hexdigest(),
        'witness_sha256': hashlib.sha256(blob).hexdigest(), 'frozen_C_exact_moment': str(fixed_value),
        'frozen_C_rejects': fixed_value < 0, 'parametric_separator': separator,
        'certified_unphysical_via_T1': fixed_value < 0 or separator is not None,
        'numerical_packages_loaded': [], 'wall_seconds': time.monotonic()-start,
        'scope': 'Exact full old-family acceptance and frozen T1 evaluation. Unphysicality is proved by this route only if a strict negative T1 expectation is accepted.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('phase', choices=('discover', 'replay'))
    args = parser.parse_args(); discover() if args.phase == 'discover' else replay()
