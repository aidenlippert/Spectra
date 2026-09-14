"""Occupation-polynomial charge-gap certificates without determinant actions.

The termwise Hermitian row envelope is below H as an operator. A Boolean
positive decomposition bounds that envelope on N_alpha=N_beta=m/2, D>=1.
LP construction is optional; exact replay imports only the standard library.
"""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path

from experiments.marginal_determinant_tree import DeterminantOracle


def add(*polys):
    out = {}
    for poly in polys:
        for mask, value in poly.items():
            out[mask] = out.get(mask, F(0)) + value
    return {m: v for m, v in out.items() if v}


def scale(poly, value):
    return {m: v * value for m, v in poly.items() if v * value}


def multiply(a, b):
    out = {}
    for x, u in a.items():
        for y, v in b.items():
            out[x | y] = out.get(x | y, F(0)) + u * v
    return {m: v for m, v in out.items() if v}


def indicator(required, occupied):
    out = {occupied: F(1)}
    empty = required ^ occupied
    while empty:
        bit = empty & -empty
        out = multiply(out, {0: F(1), bit: F(-1)})
        empty ^= bit
    return out


def envelope(data):
    oracle = DeterminantOracle(data)
    modes, particles = oracle.modes, oracle.particles
    if modes % 4 or particles != modes // 2:
        raise ValueError('Even-site half-filled sector required')
    for word in oracle.h:
        if any(sum((2*c-1) for c, i in word if i % 2 == spin) for spin in (0, 1)):
            raise ValueError('Hamiltonian must conserve each spin population')
    poly = dict(oracle.diagonal)
    for required, occupied, _, magnitude in oracle.transitions:
        poly = add(poly, scale(indicator(required, occupied), -magnitude))
    shifts = [{0: F(-particles // 2), **{1 << i: F(1) for i in range(spin, modes, 2)}}
              for spin in (0, 1)]
    charge = {0: F(-1), **{3 << i: F(1) for i in range(0, modes, 2)}}
    return oracle, poly, shifts, charge


def _mask(value, modes):
    if type(value) is not int or not 0 <= value < (1 << modes):
        raise ValueError('Invalid occupation mask')
    return value


def _atoms(items, modes, degree):
    if type(items) is not list or len(items) > 30000:
        raise ValueError('Bounded explicit indicator list required')
    out = {}
    for item in items:
        req, occ = (_mask(item[k], modes) for k in ('required', 'occupied'))
        if occ & ~req or req.bit_count() > degree:
            raise ValueError('Invalid local indicator support')
        if type(item['weight']) is not str or F(item['weight']) < 0:
            raise ValueError('Nonnegative exact indicator weight required')
        out = add(out, scale(indicator(req, occ), F(item['weight'])))
    return out


def replay(certificate):
    if certificate.get('kind') != 'valence_charge_polynomial_v1':
        raise ValueError('Unsupported charge polynomial certificate')
    degree = certificate.get('degree')
    if type(degree) is not int or degree not in (2, 3, 4):
        raise ValueError('Supported occupation degrees are 2, 3, 4')
    oracle, f, shifts, charge = envelope(certificate)
    modes = oracle.modes
    atoms = _atoms(certificate['positive_indicators'], modes, degree)
    charge_atoms = _atoms(certificate['charge_indicators'], modes, degree-2)
    ideals = certificate['number_multipliers']
    if type(ideals) is not list or len(ideals) != 2:
        raise ValueError('Two spin-number multipliers required')
    ideal = {}
    for terms, shift in zip(ideals, shifts):
        if type(terms) is not list or len(terms) > 30000:
            raise ValueError('Bounded number multiplier required')
        p = {}
        for item in terms:
            mask = _mask(item['mask'], modes)
            if mask.bit_count() > degree-1 or type(item['coefficient']) is not str:
                raise ValueError('Invalid number multiplier')
            p = add(p, {mask: F(item['coefficient'])})
        ideal = add(ideal, multiply(shift, p))
    if type(certificate['b']) is not str:
        raise ValueError('Exact bound required')
    b = F(certificate['b'])
    residual = add(f, {0: -b}, scale(atoms, -1),
                   scale(multiply(charge, charge_atoms), -1), scale(ideal, -1))
    error = sum(map(abs, residual.values()), F(0))
    return {'lower': str(b-error), 'lower_float': float(b-error),
            'residual_l1': str(error), 'residual_l1_float': float(error),
            'envelope_terms': len(f), 'residual_terms': len(residual),
            'positive_indicators': len(certificate['positive_indicators']),
            'charge_indicators': len(certificate['charge_indicators']),
            'number_multiplier_terms': sum(map(len, ideals)),
            'determinant_actions': len(oracle.cache),
            'scope': 'Exact QHQ lower bound for the complete ionic complement at half filling and equal spin populations. Termwise row-envelope operator inequality plus nonnegative local occupation indicators, D-1 localizers, number identities and exact residual norm. No determinant enumeration in construction or replay; strength is not guaranteed.'}


def propose(source, output, degree=4, denominator=10**12):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import csc_matrix
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError('Preserve previous charge certificate')
    if type(degree) is not int or degree not in (2, 3, 4) or type(denominator) is not int or denominator < 1:
        raise ValueError('Invalid proposal budget')
    data = json.loads(source.read_text())
    base = {k: data[k] for k in ('modes', 'particles', 'hamiltonian')}
    oracle, f, shifts, charge = envelope(base)
    modes = oracle.modes
    # Explicit local-feature budget, independent of determinant-sector dimension.
    if modes > 16:
        raise ValueError('LP proposer limited to sixteen modes; replay accepts larger sectors')
    masks = [sum(1 << i for i in support) for k in range(degree+1)
             for support in combinations(range(modes), k)]
    columns, labels, bounds = [{0: F(1)}], [('b',)], [(None, None)]
    for mask in masks:
        bits = [1 << i for i in range(modes) if mask & (1 << i)]
        for assignment in range(1 << len(bits)):
            occupied = sum(bit for j, bit in enumerate(bits) if assignment & (1 << j))
            poly = indicator(mask, occupied)
            columns.append(poly); labels.append(('positive', mask, occupied)); bounds.append((0, None))
            if len(bits) <= degree-2:
                columns.append(multiply(charge, poly)); labels.append(('charge', mask, occupied)); bounds.append((0, None))
        if mask.bit_count() <= degree-1:
            for spin, shift in enumerate(shifts):
                columns.append(multiply(shift, {mask: F(1)})); labels.append(('ideal', spin, mask)); bounds.append((None, None))
    rows = sorted(set(f).union(*(set(c) for c in columns)))
    index = {mask: i for i, mask in enumerate(rows)}
    ri, ci, values = [], [], []
    for j, col in enumerate(columns):
        for mask, value in col.items():
            ri.append(index[mask]); ci.append(j); values.append(float(value))
    matrix = csc_matrix((values, (ri, ci)), shape=(len(rows), len(columns)))
    objective = np.zeros(len(columns)); objective[0] = -1
    result = linprog(objective, A_eq=matrix, b_eq=[float(f.get(m, 0)) for m in rows],
                     bounds=bounds, method='highs', options={'time_limit': 120})
    if not result.success:
        raise ValueError('No optimal numerical proposal: '+result.message)
    cert = dict(base, kind='valence_charge_polynomial_v1', degree=degree,
                b=str(F(round(result.x[0]*denominator), denominator)),
                positive_indicators=[], charge_indicators=[], number_multipliers=[[], []])
    for value, label in zip(result.x[1:], labels[1:]):
        rational = F(round(value*denominator), denominator)
        if not rational:
            continue
        if label[0] == 'ideal':
            cert['number_multipliers'][label[1]].append({'mask': label[2], 'coefficient': str(rational)})
        else:
            rational = max(F(0), rational)
            key = 'positive_indicators' if label[0] == 'positive' else 'charge_indicators'
            cert[key].append({'required': label[1], 'occupied': label[2], 'weight': str(rational)})
    receipt = replay(cert)
    output.mkdir(parents=True)
    for name, obj in [('certificate', cert), ('receipt', receipt),
                      ('proposal', {'rows': len(rows), 'columns': len(columns), 'numerical_bound': float(result.x[0]),
                                    'status': result.message, 'scope': 'LP proposal; no exact optimality assertion.'})]:
        (output/(name+'.json')).write_text(json.dumps(obj, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify'); parser.add_argument('--source'); parser.add_argument('--output')
    parser.add_argument('--degree', type=int, default=4)
    args = parser.parse_args()
    print(json.dumps(replay(json.loads(Path(args.verify).read_text())) if args.verify
                     else propose(args.source, args.output, args.degree), indent=2))
