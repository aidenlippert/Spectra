"""Exact componentwise complement bounds for finite molecular Schur proofs."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_molecular_gershgorin import matrix, upper, row_bounds, prepare as scalar_prepare, replay as scalar_replay
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_enlarged_schur import solve_positive
from experiments.marginal_general_schur import moment_recurrence, resolvent_self_energy
from experiments.marginal_component_resolvent import second_order_correction, full_recurrence


def components(a, indices):
    remaining = set(indices)
    result = []
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        group, queue = [seed], [seed]
        while queue:
            i = queue.pop()
            adjacent = sorted(j for j in remaining if a[i][j] != 0)
            remaining.difference_update(adjacent)
            group.extend(adjacent)
            queue.extend(adjacent)
        result.append(sorted(group))
    return result


def strongest_edge_partition(a, indices, max_size=7):
    """Deterministic capped merges from the actual Q matrix; no quality guarantee."""
    if type(max_size) is not int or not 1 <= max_size <= 12:
        raise ValueError('Reference block cap must be an integer from 1 through 12')
    groups = {i: {i} for i in indices}
    owner = {i: i for i in indices}
    edges = sorted((-abs(a[i][j]), i, j) for i in indices for j in indices if i < j and a[i][j])
    for _, i, j in edges:
        left, right = owner[i], owner[j]
        if left == right or len(groups[left]) + len(groups[right]) > max_size:
            continue
        groups[left].update(groups.pop(right))
        for k in groups[left]:
            owner[k] = left
    return sorted((sorted(group) for group in groups.values()), key=lambda group: group[0])


def prepare(a, retained, component_lowers=None, partition=None):
    block, _, _ = scalar_prepare(a, retained)  # validates the retained budget
    q = [i for i in range(len(a)) if i not in retained]
    if partition is None:
        partition = components(a, q)
    elif (type(partition) is not list or not partition
          or any(type(group) is not list or not group or any(type(i) is not int for i in group) for group in partition)
          or sorted(i for group in partition for i in group) != q):
        raise ValueError('Reference blocks must partition every Q index exactly once')
    membership = {i: j for j, group in enumerate(partition) for i in group}
    eta = max(sum(abs(a[i][j]) for j in q if membership[i] != membership[j]) for i in q)
    groups = []
    for group in partition:
        gamma = min(row_bounds(a, group).values())
        leakage = [[sum(a[i][k] * a[k][j] for k in group) for j in retained] for i in retained]
        groups.append({'indices': group, 'lower': gamma, 'leakage': leakage,
                       'matrix': [[a[i][j] for j in group] for i in group],
                       'coupling': [[a[i][j] for j in group] for i in retained]})
    if component_lowers is not None:
        if type(component_lowers) is not list or len(component_lowers) != len(groups):
            raise ValueError('One spectral threshold per reconstructed Q component required')
        for group, value in zip(groups, component_lowers):
            lower = F(value)
            if lower < group['lower']:
                raise ValueError('Spectral threshold must not weaken Gershgorin')
            if lower > group['lower']:
                shifted = [[x - lower * (i == j) for j, x in enumerate(row)] for i, row in enumerate(group['matrix'])]
                if ldl_pivots(shifted) is None:
                    raise ValueError('Component spectral threshold failed exact positivity')
            group['lower'] = lower
    ordered_q = [i for group in partition for i in group]
    return {'block': block, 'groups': groups, 'off_block_norm': eta,
            'off_block_matrix': [[a[i][j] if membership[i] != membership[j] else F(0)
                                  for j in ordered_q] for i in ordered_q],
            'q_coupling': [[a[i][j] for j in retained] for i in ordered_q]}


def component_moments(group):
    if 'moments_data' in group:
        return group['moments_data']
    active = [i for i, row in enumerate(group['coupling']) if any(row)]
    columns = [{j: x for j, x in enumerate(group['coupling'][i]) if x} for i in active]
    def apply(current):
        result = []
        for column in current:
            output = {i: sum(x * column.get(j, F(0)) for j, x in enumerate(row))
                      for i, row in enumerate(group['matrix'])}
            result.append({i: x for i, x in output.items() if x})
        return result
    recurrence = moment_recurrence(columns, apply, max_degree=12)
    recurrence['active'] = active
    group['moments_data'] = recurrence
    return recurrence


def moment_correction(group, lower):
    recurrence = component_moments(group)
    active = recurrence['active']
    small = resolvent_self_energy(recurrence['annihilator'], recurrence['moments'], lower, size=len(active))
    result = [[F(0)] * len(group['coupling']) for _ in group['coupling']]
    for i, row in enumerate(small):
        for j, x in enumerate(row):
            result[active[i]][active[j]] = x
    return result


def schur_matrix(data, lower, exact_response=False):
    lower = F(lower)
    shift = lower + data.get('off_block_norm', F(0))
    # An uncoupled Q component can contain a lower state; never skip this gate.
    if any(group['lower'] <= shift for group in data['groups']):
        return None
    reduced = [[x - lower * (i == j) for j, x in enumerate(row)] for i, row in enumerate(data['block'])]
    if exact_response == 'second_order':
        correction = second_order_correction(data, lower)
        return [[x - y for x, y in zip(row, correction[i])] for i, row in enumerate(reduced)]
    for group in data['groups']:
        if exact_response == 'moments':
            correction = moment_correction(group, shift)
        elif exact_response:
            shifted = [[x - shift * (i == j) for j, x in enumerate(row)] for i, row in enumerate(group['matrix'])]
            solved = solve_positive(shifted, [list(row) for row in zip(*group['coupling'])])
            correction = [[sum(x * solved[k][j] for k, x in enumerate(row)) for j in range(len(reduced))]
                          for row in group['coupling']]
        else:
            correction = [[x / (group['lower'] - shift) for x in row] for row in group['leakage']]
        reduced = [[x - y for x, y in zip(row, correction[i])] for i, row in enumerate(reduced)]
    return reduced


def schur_pivots(data, lower, exact_response=False):
    reduced = schur_matrix(data, lower, exact_response)
    return None if reduced is None else ldl_pivots(reduced)


def suggest_lower(data, witness_upper, previous_lower, exact_response=False):
    import numpy as np
    a = np.array(data['block'], dtype=float)
    corrections = [(np.array(g['leakage'], dtype=float), float(g['lower'])) for g in data['groups']]
    blocks = [(np.array(g['matrix'], dtype=float), np.array(g['coupling'], dtype=float)) for g in data['groups']]
    eta = float(data.get('off_block_norm', F(0)))
    if exact_response == 'second_order':
        r = np.array(data['off_block_matrix'], dtype=float)
        w = np.array(data['q_coupling'], dtype=float)
        def action(b, rhs):
            output, offset = [], 0
            for block, _ in blocks:
                size = len(block)
                output.append(np.linalg.solve(block - b * np.eye(size), rhs[offset:offset + size]))
                offset += size
            return np.vstack(output)
    c = min(g for _, g in corrections) - eta
    def positive(b):
        if b >= c:
            return False
        reduced = a - b * np.eye(len(a))
        if exact_response == 'second_order':
            x = action(b, w)
            y = r @ x
            reduced = reduced - (w.T @ x - x.T @ y + y.T @ action(b + eta, y))
        elif exact_response:
            for block, coupling in blocks:
                reduced = reduced - coupling @ np.linalg.solve(block - (b + eta) * np.eye(len(block)), coupling.T)
        else:
            for k, gamma in corrections:
                reduced = reduced - k / (gamma - b - eta)
        return np.linalg.eigvalsh(reduced)[0] > 0
    # Each method must establish its own starting bound; no dominance is assumed.
    if schur_pivots(data, previous_lower, exact_response) is None:
        raise ValueError('Previous bound did not transfer to component correction')
    lo, hi = float(previous_lower), min(float(witness_upper), c)
    for _ in range(60):
        middle = (lo + hi) / 2
        if positive(middle): lo = middle
        else: hi = middle
    for exponent in range(10):
        lower = max(previous_lower, F(math.floor((lo - 10. ** (exponent - 10)) * 10**12), 10**12))
        if lower <= witness_upper and schur_pivots(data, lower, exact_response) is not None:
            return lower
    raise ValueError('No exact component bound accepted')


def suggest_component_lowers(data):
    import numpy as np
    result = []
    for group in data['groups']:
        estimate = np.linalg.eigvalsh(np.array(group['matrix'], dtype=float))[0]
        for exponent in range(10):
            candidate = max(group['lower'], F(math.floor((estimate - 10. ** (exponent - 10)) * 10**12), 10**12))
            shifted = [[x - candidate * (i == j) for j, x in enumerate(row)] for i, row in enumerate(group['matrix'])]
            if candidate == group['lower'] or ldl_pivots(shifted) is not None:
                result.append(rational_text(candidate))
                break
        else:
            raise ValueError('Exact component spectral export exhausted')
    return result


def replay(certificate):
    if certificate.get('kind') != 'molecular_component_schur_v1':
        raise ValueError('Unsupported component Schur certificate')
    states, a = matrix(certificate)
    chosen = certificate.get('retained_states')
    if type(chosen) is not list or any(type(s) is not int or s not in states for s in chosen):
        raise ValueError('Invalid retained bitstrings')
    bound = certificate.get('component_bound', 'gershgorin')
    response = certificate.get('response', 'scalar')
    if bound not in ('gershgorin', 'spectral') or response not in ('scalar', 'exact', 'moments', 'second_order'):
        raise ValueError('Unsupported complement method')
    thresholds = certificate.get('component_lowers')
    if (bound == 'spectral' and thresholds is None) or (bound == 'gershgorin' and thresholds is not None):
        raise ValueError('Inconsistent component threshold recipe')
    partition = certificate.get('reference_components')
    if partition is not None:
        if (type(partition) is not list or any(type(group) is not list or
                any(type(s) is not int or s not in states for s in group) for group in partition)):
            raise ValueError('Invalid reference block bitstrings')
        partition = [[states.index(s) for s in group] for group in partition]
    selection = certificate.get('reference_selection')
    if selection is not None:
        if type(selection) is not dict or set(selection) != {'method', 'max_size'} or selection['method'] != 'strongest_edges':
            raise ValueError('Unsupported automatic reference selection recipe')
        q = [i for i, state in enumerate(states) if state not in chosen]
        if partition != strongest_edge_partition(a, q, selection['max_size']):
            raise ValueError('Reference partition does not match automatic selection')
    data = prepare(a, [states.index(s) for s in chosen], thresholds, partition)
    lower = F(certificate['lower'])
    checked = schur_pivots(data, lower, response if response in ('moments', 'second_order') else response == 'exact')
    if checked is None:
        raise ValueError('Exact component Schur positivity failed')
    witness_upper = upper(a, certificate.get('independent_upper'))
    if witness_upper < lower:
        raise ValueError('Inconsistent component interval')
    result = {'lower': rational_text(lower), 'upper': rational_text(witness_upper),
            'width': rational_text(witness_upper - lower), 'width_float': float(witness_upper - lower),
            'retained_dimension': len(chosen), 'component_dimensions': [len(g['indices']) for g in data['groups']],
            'component_lowers': [rational_text(g['lower']) for g in data['groups']],
            'component_states': [[states[i] for i in g['indices']] for g in data['groups']],
            'positive_denominators': [rational_text(g['lower'] - lower - data['off_block_norm']) for g in data['groups']],
            'off_block_norm_bound': rational_text(data['off_block_norm']),
            'schur_pivots': [rational_text(x) for x in checked],
            'component_bound': bound, 'response': response,
            'scope': 'Response of reference Q blocks with certified interblock norm shift; explicit 70-state matrix; no scaling claim'}
    if response == 'moments':
        result['response_recurrences'] = [{'degree': len(component_moments(g)['annihilator']) - 1,
                                           'active_retained_columns': len(component_moments(g)['active']),
                                           'annihilator': [rational_text(x) for x in component_moments(g)['annihilator']]}
                                          for g in data['groups']]
    if response == 'second_order':
        result['response_recurrences'] = [{'degree': len(full_recurrence(g)['annihilator']) - 1,
                                           'annihilator': [rational_text(x) for x in full_recurrence(g)['annihilator']]}
                                          for g in data['groups']]
        result['scope'] = 'Second-order interblock response; rigorous norm-bounded remainder; exact small-block polynomial inverses; explicit 70-state matrix; no scaling claim'
    if selection is not None:
        result['reference_selection'] = selection
    return result


def upgrade(path, bound='gershgorin', response='scalar'):
    source = Path(path)
    if bound not in ('gershgorin', 'spectral') or response not in ('scalar', 'exact', 'moments'):
        raise ValueError('Unsupported complement method')
    suffix = '_components' + ('_spectral' if bound == 'spectral' else '') + ('_' + response if response != 'scalar' else '')
    out = source.parent.with_name(source.parent.name + suffix)
    if out.exists():
        raise ValueError('Preserve previous component export')
    started = time.monotonic()
    c = json.loads(source.read_text())
    old = scalar_replay(c)
    states, a = matrix(c)
    data = prepare(a, [states.index(s) for s in c['retained_states']])
    thresholds = suggest_component_lowers(data) if bound == 'spectral' else None
    if thresholds is not None:
        data = prepare(a, [states.index(s) for s in c['retained_states']], thresholds)
        c['component_lowers'] = thresholds
    lower = suggest_lower(data, F(old['upper']), F(old['lower']), 'moments' if response == 'moments' else response == 'exact')
    c.update(kind='molecular_component_schur_v1', lower=rational_text(lower), component_bound=bound, response=response)
    receipt = replay(c)
    receipt.update(previous_width=old['width'], previous_width_float=old['width_float'],
                   elapsed_seconds=time.monotonic() - started, source=str(source))
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(c, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('retained_dimension', 'component_dimensions', 'width_float',
                                            'previous_width_float', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upgrade')
    parser.add_argument('--verify')
    parser.add_argument('--bound', choices=('gershgorin', 'spectral'), default='gershgorin')
    parser.add_argument('--response', choices=('scalar', 'exact', 'moments'), default='scalar')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.upgrade:
        upgrade(args.upgrade, args.bound, args.response)
    else:
        parser.error('Specify --upgrade or --verify')
