"""Exact spin-changing density generators and compact SOS acceptance."""
from fractions import Fraction as F
from math import lcm
import json
import time

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import canonical, mono, product, verify
from research.joint_patterns_20260913.core import expand_certificate as baseline, pack_rows
from research.molecular_collective_20260913.core import extract, factor_operators, tail_replay
from research.spin_subspace_20260913.core import combine


def generator(patterns, ref, modes):
    if not isinstance(ref, (list, tuple)) or len(ref) != 4 or any(type(v) is not int for v in ref):
        raise ValueError('Expected [pattern, creation spin, annihilation spin, mode]')
    k, s, t, mode = ref
    if not 0 <= mode < modes:
        raise ValueError('Invalid annihilation mode')
    a = mono(((0, mode),))
    if (k, s, t) == (-1, -1, -1):
        return a
    if not 0 <= k < len(patterns) or s not in (0, 1) or t not in (0, 1):
        raise ValueError('Invalid density spin indices')
    # Extract the common spatial matrix from its up/up copy; introduce the
    # requested spin indices without assuming the component is Hermitian.
    density = {}
    for w, c in patterns[k].items():
        if len(w) != 2 or w[0][0] != 1 or w[1][0] != 0 or w[0][1] % 2 != w[1][1] % 2:
            raise ValueError('Expected spin-summed one-body spatial pattern')
        if w[0][1] % 2 == 0:
            density[((1, w[0][1]+s), (0, w[1][1]+t))] = c
    return product(density, a)


def frames(patterns, modes, signature):
    parts = {}
    for k, s, t in [(-1, -1, -1)]+[(k, s, t) for k in range(len(patterns)) for s in range(2) for t in range(2)]:
        for mode in range(modes):
            ref = [k, s, t, mode]; q = generator(patterns, ref, modes)
            keys = {signature(w) for w in q}
            if not q or len(keys) != 1:
                raise ValueError('Empty or symmetry-mixed completed spin generator')
            key = next(iter(keys))
            g = parts.setdefault(key, {'name': 'spin_complete_'+str(key), 'kind': 'anticommutator', 'generators': [], 'polynomials': []})
            g['generators'].append(ref); g['polynomials'].append(q)
    return [parts[key] for key in sorted(parts)]


def expand_certificate(data, tail, cert):
    if cert.get('kind') != 'spin_completion_subspace_v1':
        raise ValueError('Unknown completed-spin certificate')
    expanded = baseline(data, tail, {**cert, 'kind': 'joint_density_anticommutator_v1', 'anti_blocks': []})
    p = extract(data, tail['center_number']); patterns = [q for _, q in factor_operators(p, tail)]
    den = cert['denominator']; common = den; rows_by_block = []
    if not isinstance(cert['anti_blocks'], list) or len(cert['anti_blocks']) > 8:
        raise ValueError('At most eight completed-spin blocks admitted')
    for block in cert['anti_blocks']:
        refs = block['generators']; directions = block['directions']; dd = block['direction_denominator']
        if not refs or len(refs) > p['modes']*(4*len(patterns)+1):
            raise ValueError('Generator count exceeds completed-spin budget')
        polys = [generator(patterns, ref, p['modes']) for ref in refs]
        if type(dd) is not int or dd <= 0:
            raise ValueError('Positive integer direction denominator required')
        if not directions or len(directions) > len(refs) or any(len(row) != len(refs) or any(type(c) is not int for c in row) or not any(row) for row in directions):
            raise ValueError('Invalid integer subspace directions')
        factors = block['factor']
        if len(factors) > len(directions) or any(len(row) != len(directions) or any(type(c) is not int for c in row) for row in factors):
            raise ValueError('Invalid integer subspace factor')
        rows = []
        for row in factors:
            coefficients = [sum(c*direction[j] for c, direction in zip(row, directions)) for j in range(len(refs))]
            q = combine(polys, coefficients, den*dd)
            if q:
                rows.append(q); common = lcm(common, *(F(c).denominator for c in q.values()))
        if rows:
            rows_by_block.append((block['name'], rows))
    blocks = [{**b, 'factor': [[c*(common//den) for c in row] for row in b['factor']]} for b in expanded['blocks']]
    for name, rows in rows_by_block:
        blocks.append(pack_rows(name+'/B', rows, common))
        blocks.append(pack_rows(name+'/B_dagger', [canonical(adj(q)) for q in rows], common))
    return {**expanded, 'blocks': blocks, 'denominator': common}


def replay(data, tail, cert):
    start = time.monotonic(); tail_receipt = tail_replay(data, tail)
    expanded = expand_certificate(data, tail, cert); built = time.monotonic(); exact = verify(expanded)
    if exact['residual_max_degree'] > 4:
        raise AssertionError('Completed spin anticommutators left degree-six terms')
    lower = F(exact['lower'])+F(tail_receipt['lower_operator_shift_Ha'])
    return {'retained': exact, 'tail': tail_receipt, 'original_lower_Ha': str(lower),
        'original_lower_float_Ha': float(lower), 'expand_seconds': built-start,
        'SOS_replay_seconds': time.monotonic()-built, 'replay_seconds': time.monotonic()-start,
        'expanded_certificate_bytes': len(json.dumps(expanded, separators=(',', ':')).encode())+1,
        'many_body_states_enumerated': 0,
        'scope': 'Exact molecular lower with the same collective tail; no optimum claim.'}
