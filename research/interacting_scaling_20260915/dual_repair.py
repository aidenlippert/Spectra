"""Exact singlet affine/nullspace repair for an explicitly declared SOS family.

The output bounds the attainable LOWER certificate, not the physical energy.
Numerical proposals and operator-family export are separate from acceptance.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
from experiments.marginal_symbolic import add, adj, canonical, decode, encode, mono, product, scale, number_shift
from research.collective_completion_20260914.spin_replay import setup, alpha_shift
from research.collective_completion_20260914.spin_screen import spin_squared, ladder_ideal
from research.certificate_scaling.commutator_dual_witness import psd
from research.interacting_scaling_20260915.dictionary import representative
from research.interacting_scaling_20260915.singlet_trace import singlet_trace, spatial_charge
from research.sector_quotient_20260914.fast_twirl import twirl
from research.interacting_scaling_20260915.budget import dump


def export(case, proposal, output):
    import numpy as np
    from scipy import sparse
    output.mkdir(parents=True, exist_ok=False)
    prepared = case/'prepared'
    meta = json.loads((prepared/'frame.json').read_text())
    data = json.loads((case/'fixture.json').read_text())
    if meta.get('magnetization', 0):
        raise ValueError('This repair is for a singlet family')
    m, n = data['modes'], data['particles']
    bases = np.load(prepared/'bases.npz')
    groups = []
    for k, block in enumerate(meta['blocks']):
        words = meta['groups'][block['physical_group']]['words']
        V = bases[f'V_{k}']
        integers = np.rint(6*V).astype(np.int64)
        if np.max(abs(integers/6-V)) > 1e-14:
            raise ValueError('Operator basis is not the exact supported denominator-six basis')
        groups.append([encode(canonical({tuple(map(tuple, w)): F(int(integers[i, j]), 6)
            for i, w in enumerate(words) if integers[i, j]})) for j in range(V.shape[1])])
    ideals = [product(number_shift(m, n), decode(p, m, 4)) for p in meta['number_basis']]
    ideals += [product(alpha_shift(m, n), decode(p, m, 2)) for p in meta['spin_basis']]
    ideals += [spin_squared(m)]+[ladder_ideal(m, decode(p, m, 2)) for p in meta['ladder_basis']]
    h, hs, delta = setup(data)
    family = {'modes': m, 'particles': n, 'hamiltonian': encode(h), 'spin_defect_Ha': str(delta),
        'groups': groups, 'ideals': [encode(p) for p in ideals], 'spin_average': True,
        'frame_sha256': hashlib.sha256((prepared/'frame.json').read_bytes()).hexdigest(),
        'basis_file_sha256': hashlib.sha256((prepared/'bases.npz').read_bytes()).hexdigest(),
        'declared_design': json.loads((case/'design.json').read_text()),
        'scope': 'Exactly the exported polynomial-coordinate Gram family and ideals; native twirled coefficient-L1 residual.'}
    dump(output/'family.json', family)
    T = sparse.load_npz(prepared/'twirl.npz')
    selected = np.load(prepared/'selected.npy')
    z = np.load(case/proposal/'checkpoint.npz')
    values = -T[selected].T@(np.load(prepared/'scale.npy')*z['y'])/np.load(prepared/'weights.npy')
    dump(output/'proposal.json', [{'word': w, 'value': str(F(round(float(v)*10**12), 10**12))}
        for w, v in zip(meta['rows'], values)])


def compressed(poly):
    result = {}
    for w, c in twirl(poly).items():
        key = representative(w)
        result[key] = result.get(key, F(0))+c
    return {w: c for w, c in result.items() if c}


def structure(family):
    m, n = family['modes'], family['particles']
    from research.nvidia_followup_20260915.strict_replay import require_supported_sector
    require_supported_sector(family)
    if family['spin_average'] is not True:
        raise ValueError('Singlet twirl required')
    groups = [[decode(p, m, 3) for p in group] for group in family['groups']]
    ideals = [compressed(decode(p, m, 6)) for p in family['ideals']]
    grams = []
    words = {()}
    for group in groups:
        matrix = [[{} for q in group] for p in group]
        for i, p in enumerate(group):
            for j in range(i, len(group)):
                entry = compressed(product(adj(p), group[j]))
                reverse = compressed(product(adj(group[j]), p))
                if entry != reverse:
                    raise ValueError('The real symmetric Gram functional did not reproduce')
                matrix[i][j] = matrix[j][i] = entry
                words.update(entry)
        grams.append(matrix)
    objective = compressed(decode(family['hamiltonian'], m, 4))
    for p in ideals+[objective]: words.update(p)
    return groups, grams, ideals, objective, words


def evaluate(poly, moments):
    return sum((c*moments.get(w, F(0)) for w, c in poly.items()), F(0))


def eliminate(equations):
    """Exact sparse affine elimination; zero rows are checked for consistency."""
    pivots = {}
    for polynomial, rhs in equations:
        row = {w: F(v) for w, v in polynomial.items() if v}
        rhs = F(rhs)
        while row:
            key = min(row, key=lambda w: (len(w), w))
            if key not in pivots:
                value = row[key]
                pivots[key] = ({w: c/value for w, c in row.items()}, rhs/value)
                break
            old, value = pivots[key]
            coefficient = row[key]
            rhs -= coefficient*value
            for w, c in old.items():
                row[w] = row.get(w, F(0))-coefficient*c
                if not row[w]: del row[w]
        else:
            if rhs: raise ValueError('Inconsistent exact affine repair constraints')
    return pivots


def nullspace(matrix):
    n = len(matrix)
    # Reuse the exact affine elimination with integer coordinate keys encoded
    # as sortable tuples, independently of molecular coefficient identities.
    keys = [((0, j),) for j in range(n)]
    pivots = eliminate([({keys[j]: v for j, v in enumerate(row) if v}, F(0)) for row in matrix])
    vectors = []
    for free in (key for key in keys if key not in pivots):
        values = {key: F(key == free) for key in keys}
        for key, (row, rhs) in sorted(pivots.items(), key=lambda item: item[0], reverse=True):
            values[key] = rhs-sum(c*values[w] for w, c in row.items() if w != key)
        vectors.append([values[key] for key in keys])
    return vectors


def check(family, witness):
    start = time.monotonic()
    _, grams, ideals, objective, words = structure(family)
    moments = {}
    for item in witness['moments']:
        w = tuple(map(tuple, item['word']))
        if w in moments or w != representative(w) or type(item['value']) is not str:
            raise ValueError('Canonical unique exact moment required')
        moments[w] = F(item['value'])
    if set(moments) != words or moments[()] != 1 or any(abs(v) > 1 for v in moments.values()):
        raise ValueError('Moment support, normalization or coefficient box failed')
    if any(evaluate(p, moments) for p in ideals):
        raise ValueError('An exact singlet ideal equation failed')
    stats = [psd([[evaluate(p, moments) for p in row] for row in matrix]) for matrix in grams]
    value = evaluate(objective, moments)
    _, _, delta = setup(family)
    if F(family['spin_defect_Ha']) != delta:
        raise ValueError('The original-H spin allowance did not reproduce')
    return {'dual_objective_Ha': str(value), 'dual_objective_float_Ha': float(value),
        'original_H_lower_family_upper_Ha': str(value-delta),
        'Gram_checks': stats, 'ideal_equations': len(ideals), 'moment_count': len(moments),
        'replay_seconds': time.monotonic()-start,
        'scope': 'Upper bound on all b minus twirled coefficient-L1 residual attainable in the exported singlet family. Not a physical lower bound.'}


def repair(folder):
    start = time.monotonic()
    family = json.loads((folder/'family.json').read_text())
    m, n = family['modes'], family['particles']
    groups, grams, ideals, objective, words = structure(family)
    trace = {w: singlet_trace(mono(w), m, n) for w in words}
    equations = [(mono(()), F(1))]+[(p, F(0)) for p in ideals]
    null_vectors = 0
    for group, gram in zip(groups, grams):
        charges = defaultdict(list)
        for i, p in enumerate(group): charges[spatial_charge(p)].append(i)
        for positions in charges.values():
            matrix = [[evaluate(gram[i][j], trace) for j in positions] for i in positions]
            for vector in nullspace(matrix):
                null_vectors += 1
                for row in gram:
                    equation = {}
                    for j, v in zip(positions, vector):
                        if v:
                            for w, c in row[j].items(): equation[w] = equation.get(w, F(0))+v*c
                    equations.append(({w: c for w, c in equation.items() if c}, F(0)))
    for w, value in trace.items():
        if abs(value) == 1: equations.append(({w: F(1)}, value))
    if any(evaluate(p, trace) != rhs for p, rhs in equations):
        raise ValueError('Singlet trace does not satisfy the declared repair equations')
    pivots = eliminate(equations)
    proposed = {tuple(map(tuple, item['word'])): F(item['value']) for item in json.loads((folder/'proposal.json').read_text())}
    values = {w: proposed.get(w, F(0)) for w in words}
    for w, (row, rhs) in sorted(pivots.items(), key=lambda item: (len(item[0]), item[0]), reverse=True):
        values[w] = rhs-sum(c*values[v] for v, c in row.items() if v != w)
    dump(folder/'affine_repair.json', {'rank': len(pivots), 'equations': len(equations), 'forced_trace_null_vectors': null_vectors,
        'seconds': time.monotonic()-start, 'full_fixed_N_determinants_enumerated': 0})
    attempts = []
    mixtures = [F(0)]+[F(1, 10**k) for k in range(10, 0, -1)]+[F(1, 4), F(1, 2), F(3, 4), F(1)]
    for mix in mixtures:
        candidate = {w: (1-mix)*values[w]+mix*trace[w] for w in words}
        try:
            if any(abs(v) > 1 for v in candidate.values()): raise ValueError('Coefficient box failed')
            for matrix in grams: psd([[evaluate(p, candidate) for p in row] for row in matrix])
        except ValueError as error:
            attempts.append({'mixture': str(mix), 'refusal': str(error)})
            continue
        witness = {'trace_mixture': str(mix), 'moments': [{'word': w, 'value': str(candidate[w])}
            for w in sorted(words, key=lambda w: (len(w), w))]}
        dump(folder/'witness.json', witness)
        dump(folder/'construction.json', {'seconds': time.monotonic()-start, 'attempts': attempts,
            'trace_mixture': str(mix), 'numerical_objective_before_repair': float(evaluate(objective, proposed)),
            'dual_objective_after_repair': float(evaluate(objective, candidate)),
            'independent_replay_required': True})
        return
    raise ValueError('No repair accepted, including the physical singlet trace')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='action', required=True)
    e = sub.add_parser('export'); e.add_argument('case', type=Path); e.add_argument('proposal'); e.add_argument('output', type=Path)
    r = sub.add_parser('repair'); r.add_argument('folder', type=Path)
    c = sub.add_parser('check'); c.add_argument('folder', type=Path)
    a = p.parse_args()
    if a.action == 'export': export(a.case.resolve(), a.proposal, a.output.resolve())
    elif a.action == 'repair': repair(a.folder.resolve())
    else:
        result = check(json.loads((a.folder/'family.json').read_text()), json.loads((a.folder/'witness.json').read_text()))
        dump(a.folder/'replay.json', result)
        print(json.dumps(result), flush=True)
