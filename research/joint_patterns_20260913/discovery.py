"""Bounded numerical proposer for joint learned density-pattern constraints."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import time

import cvxpy as cp
import numpy as np

from experiments.marginal_symbolic import add, canonical, encode, mono, multiplier_basis, number_shift, product, scale
from research.certificate_scaling.adaptive_block_discovery import partition
from research.certificate_scaling.commutator_dictionary import export_groups
from research.certificate_scaling.direct_sparse_discovery import sparse_columns
from research.certificate_scaling.polynomial_gram_contraction import prepare, contract, whitening_transform
from research.molecular_collective_20260913.core import digest, extract, factor_operators, retained_polynomial, tail_replay
from research.joint_patterns_20260913.core import generator, prepare_anticommutators, replay


def groups_for(h, p, tail, count, coupling):
    patterns = [q for _, q in factor_operators(p, tail)]
    if type(count) is not int or not 0 <= count <= len(patterns) or coupling not in ('joint', 'separate'):
        raise ValueError('Invalid learned-pattern selection')
    base, signature, symmetry = partition(h, p['modes'], 'quadratic', True)
    groups = [{'name': g['name'], 'kind': 'square',
        'polynomials': [canonical(mono(w)) for w in g['words']]} for g in base]
    parts = {}
    selections = [range(count)] if coupling == 'joint' else [[k] for k in range(count)]
    if count:
        for selection_id, selection in enumerate(selections):
            for k in [-1, *selection]:
                for mode in range(p['modes']):
                    ref = [k, mode]; poly = generator(patterns, ref, p['modes'])
                    if not poly:
                        continue
                    sigs = {signature(w) for w in poly}
                    if len(sigs) != 1:
                        raise ValueError('Learned generator mixes exact molecular symmetries')
                    key = (selection_id, next(iter(sigs)))
                    part = parts.setdefault(key, {'name': f'anti_{coupling}_{key}',
                        'kind': 'anticommutator', 'generators': [], 'polynomials': []})
                    part['generators'].append(ref); part['polynomials'].append(poly)
    groups.extend(parts[key] for key in sorted(parts))
    return groups, signature, symmetry


def export_compact(data, tail, groups, grams, transforms, xvalues, basis, denominator=10**10):
    baseline = []; anti = []; clipped = 0.
    for g, gram, transform in zip(groups, grams, transforms):
        if g['kind'] == 'square':
            blocks, den, negative = export_groups([g], [gram], rounding=denominator)
            if denominator % den:
                raise AssertionError('Unexpected monomial baseline denominator')
            for block in blocks:
                baseline.append({**block, 'factor': [[c*(denominator//den) for c in row] for row in block['factor']]})
            clipped += negative
        else:
            ev, u = np.linalg.eigh((gram+gram.T)/2)
            clipped += float(-ev[ev < 0].sum())
            root = np.sqrt(np.maximum(ev, 0))[:, None]*u.T
            if transform is not None:
                root = root@transform
            if not np.all(np.isfinite(root)) or np.max(abs(root), initial=0)*denominator >= 2**62:
                raise ValueError('Nonfinite or excessive proposed factor')
            rows = np.rint(root*denominator).astype(np.int64)
            anti.append({'name': g['name'], 'generators': g['generators'],
                'factor': [list(map(int, row)) for row in rows if np.any(row)]})
    multiplier = add(*(scale(poly, F(round(float(v)*10**12), 10**12)) for poly, v in zip(basis, xvalues[1:])))
    return {'kind': 'joint_density_anticommutator_v1', 'fixture_sha256': digest(data),
        'tail_sha256': digest(tail), 'b': str(F(round(float(xvalues[0])*10**12), 10**12)),
        'number_multiplier': encode(multiplier), 'denominator': denominator,
        'base_blocks': baseline, 'anti_blocks': anti}, clipped


def run(fixture, tail_path, out, count, coupling='joint', budget=90., solver='CLARABEL'):
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    data = json.loads(fixture.read_text()); tail = json.loads(tail_path.read_text())
    if data['modes'] > 12 or budget <= 20 or solver not in ('CLARABEL', 'SCS'):
        raise ValueError('Invalid pass budget or solver')
    tail_replay(data, tail)
    p = extract(data, tail['center_number']); h = retained_polynomial(p, tail)
    groups, signature, symmetry = groups_for(h, p, tail, count, coupling)
    identity_sig = signature(())
    basis = [q for q in multiplier_basis(p['modes'], max_body=1) if all(signature(w) == identity_sig for w in q)]
    free_polys = [mono(())]+[product(number_shift(p['modes'], p['particles']), q) for q in basis]
    prepared = []; allwords = set(h)|{w for q in free_polys for w in q}
    preparation = []; transforms = []; conditioning = []
    for g in groups:
        fn = prepare_anticommutators if g['kind'] == 'anticommutator' else prepare
        blocks, words, stats = fn([g]); prepared.extend(blocks); allwords.update(words)
        preparation.append(stats)
        if g['kind'] == 'anticommutator':
            W, stats = whitening_transform(g['polynomials'])
            transforms.append(W); conditioning.append(stats)
        else:
            transforms.append(None); conditioning.append(None)
    if any(len(w) > 4 or sum(2*c-1 for c, _ in w) or signature(w) != identity_sig for w in allwords):
        raise AssertionError('Coefficient row escaped degree, charge, or symmetry')
    rows = sorted(allwords, key=lambda w: (len(w), w)); lookup = {w: i for i, w in enumerate(rows)}
    free = sparse_columns(free_polys, lookup)
    rhs = np.array([float(h.get(w, 0)) for w in rows])
    x = cp.Variable(len(free_polys)); rem = cp.Variable(len(rows))
    expr = free@x; grams = []; maps = []; contraction = []; indices = []
    squared = np.asarray(free.power(2).sum(axis=1)).ravel()
    for g, block, W in zip(groups, prepared, transforms):
        matrix, ix, stats = contract(block, g['polynomials'], lookup, transform=W)
        k = len(g['polynomials']); gram = cp.Variable((k, k), PSD=True)
        expr += matrix@cp.reshape(gram, (k*k,), order='C')[ix]
        grams.append(gram); maps.append(matrix); indices.append(ix); contraction.append(stats)
        squared += np.asarray(matrix.power(2).sum(axis=1)).ravel()
    row_scale = 1/np.maximum(1., np.sqrt(squared))
    equality = cp.multiply(row_scale, expr+rem-rhs) == 0
    problem = cp.Problem(cp.Minimize(cp.norm1(rem)-x[0]), [equality, rem[lookup[()]] == 0])
    built = time.monotonic()
    remaining = budget-(built-start)-20.
    if remaining < 1:
        raise RuntimeError('Construction exhausted the solve/replay reserve')
    solver_seconds = min(60., remaining)
    options = ({'tol_gap_abs': 1e-8, 'tol_gap_rel': 1e-8, 'tol_feas': 1e-8,
        'max_iter': 1000, 'time_limit': solver_seconds} if solver == 'CLARABEL' else
        {'eps': 1e-8, 'max_iters': 150000, 'time_limit_secs': solver_seconds})
    construction = {'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'pattern_constraints': count, 'retained_patterns': len(tail['factors']), 'coupling': coupling,
        'modes': p['modes'], 'particles': p['particles'], 'symmetry': symmetry,
        'base_Gram_dimensions': [len(g['polynomials']) for g in groups if g['kind'] == 'square'],
        'anti_Gram_dimensions': [len(g['polynomials']) for g in groups if g['kind'] == 'anticommutator'],
        'Gram_entries': sum(len(g['polynomials'])**2 for g in groups),
        'coefficient_rows': len(rows), 'coefficient_max_degree': max(map(len, rows)),
        'gram_map_nonzeros': sum(matrix.nnz for matrix in maps), 'ideal_variables': len(basis),
        'preparation': preparation, 'contraction': contraction, 'conditioning': conditioning,
        'construction_seconds': built-start, 'solver': solver, 'solver_options': options,
        'budget_seconds': budget, 'replay_reserve_seconds': 20,
        'many_body_states_enumerated': 0}
    (out/'construction.json').write_text(json.dumps(construction, indent=2)+'\n')
    print(json.dumps({'stage': 'solve', 'patterns': count, 'coupling': coupling,
        'anti_dimensions': construction['anti_Gram_dimensions'], 'construction_seconds': built-start}), flush=True)
    problem.solve(solver=solver, **options)
    solved = time.monotonic()
    if x.value is None or any(q.value is None for q in grams):
        raise RuntimeError('Solver did not return exportable factors')
    values = [q.value for q in grams]
    np.savez_compressed(out/'proposal.npz', x=x.value,
        **{f'gram_{i}': v for i, v in enumerate(values)},
        **{f'transform_{i}': W for i, W in enumerate(transforms) if W is not None})
    if equality.dual_value is not None:
        (out/'dual_proposal.json').write_text(json.dumps({'rows': rows,
            'values': [float(v) for v in equality.dual_value*row_scale],
            'scope': 'Floating proposal; not an accepted cone obstruction.'}, separators=(',', ':'))+'\n')
    cert, clipped = export_compact(data, tail, groups, values, transforms, x.value, basis)
    encoded = json.dumps(cert, separators=(',', ':'))+'\n'
    exported = time.monotonic()
    accepted = replay(data, tail, cert)
    verified = time.monotonic()
    (out/'certificate.json').write_text(encoded)
    result = {**construction, 'status': problem.status,
        'floating_lower_Ha': -float(problem.value), 'clipped_negative_Gram_mass': clipped,
        'solve_seconds': solved-built, 'export_seconds': exported-solved,
        'accept_seconds': verified-exported, 'internal_wall_seconds': time.monotonic()-start,
        'compact_certificate_bytes': len(encoded.encode()),
        'compact_factor_nonzeros': sum(sum(c != 0 for c in row) for b in cert['base_blocks']+cert['anti_blocks'] for row in b['factor']),
        'accepted': accepted}
    (out/'receipt.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'stage': 'accepted', 'original_lower_Ha': accepted['original_lower_float_Ha'],
        'wall_seconds': result['internal_wall_seconds'], 'status': problem.status}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fixture', type=Path, required=True); parser.add_argument('--tail', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); parser.add_argument('--count', type=int, required=True)
    parser.add_argument('--coupling', choices=['joint', 'separate'], default='joint')
    parser.add_argument('--budget', type=float, default=90.)
    parser.add_argument('--solver', choices=['CLARABEL', 'SCS'], default='CLARABEL')
    args = parser.parse_args()
    run(args.fixture, args.tail, args.out, args.count, args.coupling, args.budget, args.solver)
