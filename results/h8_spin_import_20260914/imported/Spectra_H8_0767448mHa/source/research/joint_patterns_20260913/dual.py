"""Exact feasible dual ceiling for the ten-pattern H6 anticommutator cone.

No numerical library is imported. The saved floating proposal is untrusted;
only affine identities, coefficient bounds, and exact PSD tests accept it.
"""
from fractions import Fraction as F
from itertools import combinations
from math import lcm
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import canonical, mono, multiplier_basis, number_shift, product, validate_word, word_product
from research.certificate_scaling.commutator_dual_witness import evaluate, groups as quadratic_groups, moment_decode, psd, seed
from research.molecular_collective_20260913.core import digest, extract, factor_operators, retained_polynomial, tail_replay
from research.joint_patterns_20260913.core import generator, prepare_anticommutators


def cone_groups(p, tail, masks):
    m = p['modes']; h = retained_polynomial(p, tail)
    if m != 12 or p['particles'] != 6 or len(tail['factors']) != 10:
        raise ValueError('This obstruction targets the frozen ten-pattern H6 cone')
    if not isinstance(masks, list) or len(masks) > m or any(type(mask) is not int or not 0 < mask < 1 << m for mask in masks):
        raise ValueError('Invalid parity symmetry masks')
    if any(sum(2*c-1 for c, i in w if i % 2 == 0) for w in h):
        raise ValueError('Alpha-spin charge is not a symmetry')
    if any(sum((mask >> i) & 1 for _, i in w) % 2 for mask in masks for w in h):
        raise ValueError('Claimed parity does not preserve the retained Hamiltonian')
    def signature(w):
        return (sum(2*c-1 for c, i in w if i % 2 == 0),)+tuple(sum((mask >> i) & 1 for _, i in w) % 2 for mask in masks)
    def split(polys, kind):
        parts = {}
        for poly in polys:
            keys = {signature(w) for w in poly}
            if len(keys) != 1:
                raise ValueError('Generator does not have definite symmetry')
            parts.setdefault(next(iter(keys)), []).append(poly)
        return [{'kind': kind, 'polynomials': parts[key]} for key in sorted(parts)]
    result = []
    for polys in quadratic_groups(h, m, False, False):
        result.extend(split(polys, 'square'))
    patterns = [q for _, q in factor_operators(p, tail)]
    polys = [generator(patterns, [k, mode], m) for k in range(-1, 10) for mode in range(m)]
    result.extend(split(polys, 'anticommutator'))
    return result, signature


def integer_grams(polys, moments_list, kind):
    """Integer congruence avoids Fraction work in dense polynomial expansion.

    All returned matrices share the same positive scale: coefficient_den^2
    times moment_den. Thus they can be mixed before an exact PSD check.
    """
    words = sorted({w for p in polys for w in p}, key=lambda w: (len(w), w))
    lookup = {w: i for i, w in enumerate(words)}; size = len(words)
    coefficient_den = lcm(*(F(c).denominator for p in polys for c in p.values()))
    moment_den = lcm(*(v.denominator for y in moments_list for v in y.values()))
    columns = [[(lookup[w], int(c*coefficient_den)) for w, c in p.items()] for p in polys]
    values = [{w: int(v*moment_den) for w, v in y.items()} for y in moments_list]
    matrices = [[0]*(size*size) for _ in moments_list]
    if kind == 'anticommutator':
        prepared, _, _ = prepare_anticommutators([{'polynomials': polys}])
        terms = prepared[0]['terms']
    elif kind == 'square':
        terms = []
        for i, left in enumerate(words):
            dagger = tuple((1-c, j) for c, j in reversed(left))
            for j, right in enumerate(words):
                for w, c in word_product(dagger, right):
                    if len(w) > 4:
                        raise ValueError('Quadratic cone exceeded its degree budget')
                    terms.append((w, i*size+j, c))
    else:
        raise ValueError('Unknown dual Gram kind')
    for w, index, c in terms:
        for matrix, y in zip(matrices, values):
            matrix[index] += c*y.get(w, 0)
    result = []
    for matrix in matrices:
        if any(matrix[i*size+j] != matrix[j*size+i] for i in range(size) for j in range(i)):
            raise ValueError('Non-Hermitian monomial moment matrix')
        left = [[sum(c*matrix[w*size+j] for w, c in column) for j in range(size)] for column in columns]
        gram = [[sum(left[i][w]*c for w, c in column) for column in columns] for i in range(len(polys))]
        result.append(gram)
    return result, coefficient_den**2*moment_den


def affine_round(raw, m, n):
    if len(raw['rows']) != len(raw['values']):
        raise ValueError('Dual row/value count mismatch')
    proposed = {validate_word(w, m, 4): F(round(v*10**12), 10**12) for w, v in zip(raw['rows'], raw['values'])}
    if len(proposed) != len(raw['rows']):
        raise ValueError('Repeated proposed moment')
    representatives = {}
    for w in proposed:
        if canonical(mono(w)) != mono(w) or sum(2*c-1 for c, _ in w):
            raise ValueError('Noncanonical or unbalanced proposed moment')
        partner = canonical(adj(mono(w)))
        if len(partner) != 1 or next(iter(partner.values())) != 1:
            raise ValueError('Unexpected balanced adjoint')
        representatives[w] = min(w, next(iter(partner)))
    keys = sorted(set(representatives.values()), key=lambda w: (len(w), w)); lookup = {w: i for i, w in enumerate(keys)}
    values = [(proposed.get(w, F(0))+proposed.get(next(iter(canonical(adj(mono(w))))), F(0)))/2 for w in keys]
    ideals = [mono(())]+[product(number_shift(m, n), p) for p in multiplier_basis(m, max_body=1)]
    pivots = {}
    for index, poly in enumerate(ideals):
        row = {}; rhs = F(index == 0)
        for w, c in poly.items():
            if w in representatives:
                j = lookup[representatives[w]]; row[j] = row.get(j, F(0))+c
        row = {j: v for j, v in row.items() if v}
        while row:
            j = min(row); v = row[j]
            if j not in pivots:
                pivots[j] = ({k: c/v for k, c in row.items()}, rhs/v)
                break
            old, b = pivots[j]; rhs -= v*b
            for k, c in old.items():
                row[k] = row.get(k, F(0))-v*c
                if not row[k]:
                    del row[k]
        else:
            if rhs:
                raise ValueError('Inconsistent affine dual constraints')
    for j, (row, rhs) in sorted(pivots.items(), reverse=True):
        values[j] = rhs-sum(c*values[k] for k, c in row.items() if k != j)
    return {w: values[lookup[r]] for w, r in representatives.items()}, len(pivots)


def check(data, tail, witness):
    start = time.monotonic()
    if witness.get('kind') != 'joint_density_dual_v1' or witness.get('fixture_sha256') != digest(data) or witness.get('tail_sha256') != digest(tail):
        raise ValueError('Dual molecule/pattern binding failed')
    tail_result = tail_replay(data, tail)
    p = extract(data, tail['center_number']); groups, signature = cone_groups(p, tail, witness['parity_masks'])
    y = moment_decode(witness['moments'], p['modes']); zero = signature(())
    if y.get(()) != 1 or any(abs(v) > 1 for v in y.values()):
        raise ValueError('Dual normalization or coefficient box failed')
    if any(len(w) > 4 or sum(2*c-1 for c, _ in w) or (v and signature(w) != zero) for w, v in y.items()):
        raise ValueError('Dual moment degree, charge, or symmetry mismatch')
    if any(evaluate(canonical(adj(mono(w))), y) != v for w, v in y.items()):
        raise ValueError('Dual functional is not Hermitian')
    ideals = [product(number_shift(p['modes'], p['particles']), q) for q in multiplier_basis(p['modes'], max_body=1)]
    if any(evaluate(q, y) for q in ideals):
        raise ValueError('Dual number-ideal equality failed')
    stats = []
    for group in groups:
        matrices, _ = integer_grams(group['polynomials'], [y], group['kind'])
        stats.append({'kind': group['kind'], **psd(matrices[0])})
    ceiling = evaluate(retained_polynomial(p, tail), y)
    original_ceiling = ceiling+F(tail_result['lower_operator_shift_Ha'])
    return {'retained_lower_ceiling_Ha': str(ceiling), 'retained_lower_ceiling_float_Ha': float(ceiling),
        'original_lower_ceiling_Ha': str(original_ceiling), 'original_lower_ceiling_float_Ha': float(original_ceiling),
        'exact_PSD_checks': stats, 'number_ideal_equalities': len(ideals), 'moments': len(y),
        'replay_seconds': time.monotonic()-start, 'many_body_states_enumerated': 0,
        'scope': 'Ceiling on every b minus coefficient-l1 residual lower bound in the full quadratic plus joint ten-density-pattern anticommutator cone, with the body-one number ideal and fixed lower tail shift. This is not a physical energy lower bound.'}


def propose(case, fixture, tail_path, out):
    start = time.monotonic(); out.mkdir(parents=True, exist_ok=False)
    data = json.loads(fixture.read_text()); tail = json.loads(tail_path.read_text())
    raw_path = case/'dual_proposal.json'; raw = json.loads(raw_path.read_text())
    construction = json.loads((case/'construction.json').read_text())
    if construction['pattern_constraints'] != 10 or construction['coupling'] != 'joint' or construction['fixture_sha256'] != digest(data) or construction['tail_sha256'] != digest(tail):
        raise ValueError('Wrong source cone for this diagnostic')
    p = extract(data, tail['center_number']); masks = construction['symmetry']['parity_masks']
    groups, _ = cone_groups(p, tail, masks)
    expected_base = construction['base_Gram_dimensions']; expected_anti = construction['anti_Gram_dimensions']
    if [len(g['polynomials']) for g in groups if g['kind'] == 'square'] != expected_base or [len(g['polynomials']) for g in groups if g['kind'] == 'anticommutator'] != expected_anti:
        raise ValueError('Rebuilt exact cone dimensions differ from discovery')
    rounded, rank = affine_round(raw, p['modes'], p['particles']); allwords = set(rounded)
    for k in range(3):
        for inds in combinations(range(p['modes']), k):
            allwords.add(tuple((1, i) for i in inds)+tuple((0, i) for i in inds))
    trace = {w: seed(w, p['modes'], p['particles']) for w in allwords}
    pairs = [integer_grams(g['polynomials'], [rounded, trace], g['kind'])[0] for g in groups]
    attempts = []
    for mix in (F(0), F(1, 10**8), F(1, 10**7), F(1, 10**6), F(1, 10**5), F(1, 10**4), F(1, 1000)):
        try:
            for candidate, seed_matrix in pairs:
                matrix = [[(mix.denominator-mix.numerator)*x+mix.numerator*z for x, z in zip(row, seed_row)] for row, seed_row in zip(candidate, seed_matrix)]
                psd(matrix)
        except ValueError as error:
            attempts.append({'mix': str(mix), 'refusal': str(error)})
            continue
        y = {w: (1-mix)*rounded.get(w, F(0))+mix*trace[w] for w in allwords}
        witness = {'kind': 'joint_density_dual_v1', 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
            'parity_masks': masks, 'trace_mixture': str(mix),
            'moments': [{'word': w, 'value': str(v)} for w, v in sorted(y.items(), key=lambda item: (len(item[0]), item[0])) if v or not w]}
        receipt = check(data, tail, witness)
        receipt.update({'trace_mixture': str(mix), 'affine_rank': rank, 'attempts': attempts,
            'source_proposal_sha256': hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            'wall_seconds': time.monotonic()-start})
        (out/'witness.json').write_text(json.dumps(witness, separators=(',', ':'))+'\n')
        (out/'receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
        print(json.dumps(receipt), flush=True); return receipt
    (out/'refusals.json').write_text(json.dumps(attempts, indent=2)+'\n')
    raise ValueError('No exact dual repair accepted within the frozen mixture list')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--case', type=Path)
    parser.add_argument('--fixture', type=Path, required=True); parser.add_argument('--tail', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True); parser.add_argument('--witness', type=Path)
    args = parser.parse_args()
    if bool(args.case) == bool(args.witness):
        parser.error('Choose proposal or fresh witness replay')
    if args.witness:
        result = check(json.loads(args.fixture.read_text()), json.loads(args.tail.read_text()), json.loads(args.witness.read_text()))
        if any(name in sys.modules for name in ('numpy', 'scipy', 'cvxpy', 'pyscf')):
            raise AssertionError('Numerical package loaded on dual acceptance path')
        args.out.write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result), flush=True)
    else:
        propose(args.case, args.fixture, args.tail, args.out)
