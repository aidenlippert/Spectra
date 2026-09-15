"""Exact compact certificates for subspaces of the learned spin generators."""
from fractions import Fraction as F
from math import lcm
import json
import time

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import add, canonical, scale, verify
from research.joint_patterns_20260913.core import expand_certificate as expand_baseline, pack_rows
from research.joint_patterns_20260913.spin_diagnostic import spin_generator
from research.molecular_collective_20260913.core import extract, factor_operators, tail_replay


def combine(polys, coefficients, denominator):
    return add(*(scale(p, F(c, denominator)) for p, c in zip(polys, coefficients) if c))


def expand_certificate(data, tail, cert):
    if cert.get('kind') != 'spin_pattern_subspace_v1':
        raise ValueError('Unknown spin-subspace certificate')
    # Reuse every old binding, sector, baseline degree and multiplier gate.
    expanded = expand_baseline(data, tail, {**cert, 'kind': 'joint_density_anticommutator_v1', 'anti_blocks': []})
    p = extract(data, tail['center_number']); patterns = [q for _, q in factor_operators(p, tail)]
    den = cert['denominator']; common = den; rows_by_block = []
    if not isinstance(cert['anti_blocks'], list) or len(cert['anti_blocks']) > 4:
        raise ValueError('At most four spin-subspace blocks admitted')
    for block in cert['anti_blocks']:
        refs = block['generators']; directions = block['directions']; dd = block['direction_denominator']
        if not refs or len(refs) > p['modes']*(2*len(patterns)+1):
            raise ValueError('Generator count exceeds the spin-frame budget')
        polys = [spin_generator(patterns, ref, p['modes']) for ref in refs]
        if type(dd) is not int or dd <= 0:
            raise ValueError('Positive integer direction denominator required')
        if not directions or len(directions) > len(refs) or any(len(row) != len(refs) or any(type(c) is not int for c in row) or not any(row) for row in directions):
            raise ValueError('Invalid integer subspace directions')
        factors = block['factor']
        if len(factors) > len(directions) or any(len(row) != len(directions) or any(type(c) is not int for c in row) for row in factors):
            raise ValueError('Invalid integer subspace factor')
        rows = []
        for row in factors:
            # Compose the two small coordinate maps as integers before CAR
            # expansion. Its adjoint uses this exact same rounded row.
            coefficients = [sum(c*direction[j] for c, direction in zip(row, directions)) for j in range(len(refs))]
            q = combine(polys, coefficients, den*dd)
            if q:
                rows.append(q)
                common = lcm(common, *(F(c).denominator for c in q.values()))
        if rows:
            rows_by_block.append((block['name'], rows))
    blocks = [{**b, 'factor': [[c*(common//den) for c in row] for row in b['factor']]} for b in expanded['blocks']]
    for name, rows in rows_by_block:
        blocks.append(pack_rows(name+'/B', rows, common))
        blocks.append(pack_rows(name+'/B_dagger', [canonical(adj(q)) for q in rows], common))
    return {**expanded, 'blocks': blocks, 'denominator': common}


def replay(data, tail, cert):
    start = time.monotonic(); tail_receipt = tail_replay(data, tail)
    expanded = expand_certificate(data, tail, cert); built = time.monotonic()
    exact = verify(expanded)
    if exact['residual_max_degree'] > 4:
        raise AssertionError('Tied anticommutators left degree-six terms')
    lower = F(exact['lower'])+F(tail_receipt['lower_operator_shift_Ha'])
    return {'retained': exact, 'tail': tail_receipt,
        'original_lower_Ha': str(lower), 'original_lower_float_Ha': float(lower),
        'expand_seconds': built-start, 'SOS_replay_seconds': time.monotonic()-built,
        'replay_seconds': time.monotonic()-start,
        'expanded_certificate_bytes': len(json.dumps(expanded, separators=(',', ':')).encode())+1,
        'many_body_states_enumerated': 0,
        'scope': 'Exact lower bound for the frozen molecule. No optimality or compact upper-discovery claim.'}
