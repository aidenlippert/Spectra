"""Residual-certified directional response from sparse molecular actions."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_general_schur import gram
from experiments.marginal_enlarged_schur import solve_positive
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_sparse_molecular import replay as scalar_replay


def q_action(oracle, retained, vector):
    result = {}
    for state, amplitude in vector.items():
        if amplitude:
            for destination, coefficient in oracle.action(state).items():
                if destination not in retained:
                    result[destination] = result.get(destination, 0) + amplitude * coefficient
    return {s: x for s, x in result.items() if x}


def parse_basis(oracle, retained, recipe, limit=16):
    if limit not in (16, 32) or type(recipe) is not list or not 1 <= len(recipe) <= limit:
        raise ValueError('Response direction count exceeds its fixed budget')
    result = []
    for column in recipe:
        if type(column) is not dict:
            raise ValueError('Invalid sparse response direction')
        states, amplitudes = column.get('states'), column.get('amplitudes')
        if (type(states) is not list or not 1 <= len(states) <= 4096
                or any(not oracle.valid_state(s) or s in retained for s in states)
                or len(set(states)) != len(states) or type(amplitudes) is not list
                or len(amplitudes) != len(states) or any(type(x) is not int for x in amplitudes)
                or not any(amplitudes)):
            raise ValueError('Response directions must be nonzero integer vectors entirely in Q')
        result.append({s: F(x) for s, x in zip(states, amplitudes) if x})
    return result


def prepare(oracle, p, gamma, basis):
    retained = oracle.retained(p)
    block, leakage, _ = oracle.retained_data(p, gamma)
    w = [{s: x for s, x in oracle.action(state).items() if s not in retained} for state in p]
    images = [q_action(oracle, retained, column) for column in basis]
    metric = gram(basis, basis)
    if basis and ldl_pivots(metric) is None:
        raise ValueError('Response basis is linearly dependent')
    return {'block': block, 'leakage': leakage, 'gamma': F(gamma), 'w': w,
            'metric': metric, 'hamiltonian': gram(basis, images), 'squared': gram(images, images),
            'coupling': gram(basis, w), 'h_coupling': gram(images, w)}


def response_matrix(data, lower, numerical=False):
    gamma = float(data['gamma']) if numerical else data['gamma']
    delta = gamma - lower
    if delta <= 0:
        return None
    count = len(data['block'])
    if numerical:
        import numpy as np
        block = np.array(data['block'], dtype=float)
        leakage = np.array(data['leakage'], dtype=float)
        correction = np.zeros_like(leakage)
        if data['metric']:
            d = (np.array(data['squared'], dtype=float) - (lower + gamma) * np.array(data['hamiltonian'], dtype=float)
                 + lower * gamma * np.array(data['metric'], dtype=float))
            l = np.array(data['h_coupling'], dtype=float) - gamma * np.array(data['coupling'], dtype=float)
            correction = l.T @ np.linalg.solve(d, l)
        return block - lower * np.eye(count) - (leakage - correction) / delta
    correction = [[F(0)] * count for _ in range(count)]
    if data['metric']:
        d = [[x - (lower + gamma) * data['hamiltonian'][i][j] + lower * gamma * data['metric'][i][j]
              for j, x in enumerate(row)] for i, row in enumerate(data['squared'])]
        l = [[x - gamma * data['coupling'][i][j] for j, x in enumerate(row)] for i, row in enumerate(data['h_coupling'])]
        solved = solve_positive(d, l)
        correction = [[sum(l[k][i] * solved[k][j] for k in range(len(l))) for j in range(count)] for i in range(count)]
    return [[x - lower * (i == j) - (data['leakage'][i][j] - correction[i][j]) / delta
             for j, x in enumerate(row)] for i, row in enumerate(data['block'])]


def replay_rank_obstruction(certificate):
    if certificate.get('kind') != 'response_rank_obstruction_v1':
        raise ValueError('Unsupported response rank obstruction')
    oracle = DeterminantOracle(certificate)
    p = certificate.get('retained_states')
    oracle.retained(p)
    gamma, lower = F(certificate['complement_lower']), F(certificate['target_lower'])
    indices = certificate.get('negative_principal_indices')
    if (gamma <= lower or type(indices) is not list or not indices
            or any(type(i) is not int or not 0 <= i < len(p) for i in indices)
            or len(set(indices)) != len(indices)):
        raise ValueError('Positive denominator and distinct retained indices required')
    scalar = response_matrix(prepare(oracle, p, gamma, []), lower)
    pivots = ldl_pivots([[-scalar[i][j] for j in indices] for i in indices])
    if pivots is None:
        raise ValueError('Claimed negative principal subspace is not negative definite')
    return {'target_lower': rational_text(lower), 'complement_lower': rational_text(gamma),
            'minimum_response_dimension': len(indices),
            'negative_principal_dimension': len(indices),
            'unique_action_states': len(oracle.cache), 'referenced_determinants': oracle.referenced_state_count(),
            'scope': 'At this fixed target and complement threshold, any PSD improvement of the scalar Schur bound must have rank at least this dimension; hence the optimized residual-response family needs at least this many directions. Not a universal certificate-size lower bound or physical gap proof.'}


def replay(certificate, *, audit=None):
    kind = certificate.get('kind')
    if kind not in ('sparse_response_schur_v1', 'factor_response_schur_v1'):
        raise ValueError('Unsupported sparse response certificate')
    oracle = DeterminantOracle(certificate)
    p = certificate.get('retained_states')
    retained = oracle.retained(p)
    gamma, lower = F(certificate['complement_lower']), F(certificate['lower'])
    if gamma <= lower:
        raise ValueError('Positive complement gap required')
    if kind == 'sparse_response_schur_v1':
        if type(certificate.get('complement_tree')) is not list:
            raise ValueError('Explicit coverage tree required')
        _, coverage = oracle.cover(p, gamma, certificate['complement_tree'])
        scope = 'Optimized residual upper response with exact sparse actions and complete occupation-tree gap proof; no full Q inverse or sector matrix; no compact-frontier or scaling guarantee'
    else:
        from experiments.marginal_h6_complement import verify_complement
        coverage = verify_complement(oracle, p, gamma, certificate.get('blocks'))
        scope = 'Exact global interval from residual-response Schur positivity and complete explicit Q block factor proof; every complement coordinate is represented; no compact-frontier or scaling guarantee'
    basis = parse_basis(oracle, retained, certificate.get('response_basis'), 32 if kind == 'factor_response_schur_v1' else 16)
    data = prepare(oracle, p, gamma, basis)
    pivots = ldl_pivots(response_matrix(data, lower))
    if pivots is None:
        raise ValueError('Exact residual-response Schur positivity failed')
    upper = oracle.upper(certificate.get('independent_upper'))
    if lower > upper:
        raise ValueError('Inconsistent sparse response interval')
    if audit is not None:
        audit['sources'].update(oracle.cache)
        audit['referenced'].update(oracle.cache)
        audit['referenced'].update(s for row in oracle.cache.values() for s in row)
    return {'lower': rational_text(lower), 'upper': rational_text(upper),
            'width': rational_text(upper - lower), 'width_float': float(upper - lower),
            'retained_dimension': len(p), 'response_dimension': len(basis),
            'response_support_sizes': [len(column) for column in basis],
            'positive_denominator': rational_text(gamma - lower), 'complement_coverage': coverage,
            'unique_action_states': len(oracle.cache), 'referenced_determinants': oracle.referenced_state_count(),
            'schur_pivots': [rational_text(x) for x in pivots],
            'scope': scope}


def approximate_response(oracle, retained, rhs, lower):
    """Numerical CG proposal using only Q action; replay never invokes this."""
    def dot(left, right):
        return sum(x * right.get(s, 0.) for s, x in left.items())
    def combine(left, right, factor):
        result = dict(left)
        for s, x in right.items():
            result[s] = result.get(s, 0.) + factor * x
        return {s: x for s, x in result.items() if x}
    def action(vector):
        image = {s: float(x) for s, x in q_action(oracle, retained, vector).items()}
        return combine(image, vector, -lower)
    x, residual, direction = {}, dict(rhs), dict(rhs)
    norm = dot(residual, residual)
    initial = norm
    if not initial:
        raise ValueError('No response in the proposed retained direction')
    for step in range(256):
        image = action(direction)
        denominator = dot(direction, image)
        if denominator <= 0:
            raise ValueError('Numerical response proposal lost positive curvature')
        alpha = norm / denominator
        x = combine(x, direction, alpha)
        residual = combine(residual, image, -alpha)
        next_norm = dot(residual, residual)
        if next_norm <= initial * 1e-24:
            return x, step + 1
        direction = combine(residual, direction, next_norm / norm)
        norm = next_norm
    raise ValueError('Numerical response iteration budget exhausted')


def propose_direction(oracle, retained, w, coefficients, lower, basis):
    """Produce one integer response direction; callers must replay its proof."""
    rhs = {}
    for coefficient, column in zip(coefficients, w):
        for s, value in column.items():
            rhs[s] = rhs.get(s, 0.) + float(value) * coefficient
    column, steps = approximate_response(oracle, retained, rhs, lower)
    for old_column in basis:
        norm = sum(float(x) ** 2 for x in old_column.values())
        overlap = sum(float(x) * column.get(s, 0.) for s, x in old_column.items()) / norm
        for s, x in old_column.items():
            column[s] = column.get(s, 0.) - overlap * float(x)
    norm = math.sqrt(sum(x*x for x in column.values()))
    if not norm:
        raise ValueError('Response direction did not enlarge the span')
    integers = {s: int(round(x / norm * 10**10)) for s, x in column.items()}
    states = sorted(s for s, x in integers.items() if x)
    return {'states': states, 'amplitudes': [integers[s] for s in states]}, steps


def refine(path, target=F(1, 10**7)):
    source = Path(path)
    out = source.parent.with_name(source.parent.name + '_response')
    if out.exists():
        raise ValueError('Preserve previous sparse response export')
    c = json.loads(source.read_text())
    old = scalar_replay(c)
    c['kind'] = 'sparse_response_schur_v1'
    return refine_response(c, out, old, target, str(source))


def refine_factor(path, upper_path, target=F(1, 10**7)):
    from experiments.marginal_h6_complement import replay_factor
    source = Path(path)
    out = source.parent.with_name(source.parent.name + '_response')
    c = json.loads(source.read_text())
    replay_factor(c)
    upper_source = json.loads(Path(upper_path).read_text())
    # Only import a variational vector; recompute its energy in this Hamiltonian.
    c['independent_upper'] = upper_source['independent_upper']
    oracle = DeterminantOracle(c)
    upper = oracle.upper(c['independent_upper'])
    gamma = F(c['complement_lower'])
    data = prepare(oracle, c['retained_states'], gamma, [])
    for exponent in range(32):
        lower = F(math.floor(min(upper, gamma) - 2**exponent))
        if ldl_pivots(response_matrix(data, lower)) is not None:
            break
    else:
        raise ValueError('No exact initial scalar Schur bound')
    c.update(kind='factor_response_schur_v1', lower=rational_text(lower))
    old = {'upper': rational_text(upper), 'width_float': float(upper - lower)}
    return refine_response(c, out, old, target, str(source))


def refine_response(c, out, old, target, source, *, oracle=None, check=None):
    import numpy as np
    if out.exists():
        raise ValueError('Preserve previous response export')
    oracle = DeterminantOracle(c) if oracle is None else oracle
    check = replay if check is None else check
    p = c['retained_states']
    retained, gamma, upper = oracle.retained(p), F(c['complement_lower']), F(old['upper'])
    budget = 32 if c['kind'] in ('factor_response_schur_v1', 'spin_zero_response_recipe_v1') else 16
    recipe, history = list(c.get('response_basis', [])), []
    basis = parse_basis(oracle, retained, recipe, budget) if recipe else []
    started = time.monotonic()
    desired = float(upper) - 1e-10
    previous = F(c['lower'])
    best = None
    audit = {'sources': set(), 'referenced': set()}
    for iteration in range(len(basis), budget + 1):
        data = prepare(oracle, p, gamma, basis)
        if basis:
            lo, hi = float(previous), float(upper)
            for _ in range(60):
                middle = (lo + hi) / 2
                if np.linalg.eigvalsh(response_matrix(data, middle, True))[0] > 0:
                    lo = middle
                else:
                    hi = middle
            for margin in range(10):
                lower = max(previous, F(math.floor((lo - 10. ** (margin - 10)) * 10**12), 10**12))
                if ldl_pivots(response_matrix(data, lower)) is not None:
                    break
            else:
                raise ValueError('Exact residual-response export failed')
            c.update(response_basis=list(recipe), lower=rational_text(lower))
            checked = check(c, audit=audit)
            best = dict(c), checked
            previous = lower
            history.append({'response_dimension': len(basis), 'width_float': checked['width_float']})
            print(json.dumps(history[-1]), flush=True)
            if upper - lower <= target:
                break
        if iteration == budget:
            break
        _, vectors = np.linalg.eigh(response_matrix(data, desired, True))
        column, steps = propose_direction(oracle, retained, data['w'], vectors[:, 0], desired, basis)
        recipe.append(column)
        basis = parse_basis(oracle, retained, recipe, budget)
        if history:
            history[-1]['next_proposal_cg_steps'] = steps
    if best is None:
        raise ValueError('No residual-response certificate accepted')
    certificate, receipt = best
    audit['sources'].update(oracle.cache)
    audit['referenced'].update(oracle.cache)
    audit['referenced'].update(s for row in oracle.cache.values() for s in row)
    receipt.update(previous_width_float=old['width_float'], source=str(source),
                   elapsed_seconds=time.monotonic() - started,
                   proposal_unique_action_states=len(oracle.cache),
                   proposal_referenced_determinants=oracle.referenced_state_count(),
                   construction_unique_action_states=len(audit['sources']),
                   construction_referenced_determinants=len(audit['referenced']),
                   stopping_reason='target_width' if F(receipt['width']) <= target else 'response_dimension_budget')
    out.mkdir(parents=True)
    for filename, value in (('certificate.json', certificate), ('receipt.json', receipt), ('history.json', history)):
        (out / filename).write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('width_float', 'response_dimension', 'elapsed_seconds', 'stopping_reason')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refine')
    parser.add_argument('--verify')
    parser.add_argument('--refine-factor')
    parser.add_argument('--upper')
    parser.add_argument('--continue-response')
    parser.add_argument('--verify-rank-obstruction')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.refine:
        refine(args.refine)
    elif args.refine_factor and args.upper:
        refine_factor(args.refine_factor, args.upper)
    elif args.continue_response:
        source = Path(args.continue_response)
        c = json.loads(source.read_text())
        refine_response(c, source.parent.with_name(source.parent.name + '_continued'), replay(c), F(1, 10**7), str(source))
    elif args.verify_rank_obstruction:
        print(json.dumps(replay_rank_obstruction(json.loads(Path(args.verify_rank_obstruction).read_text())), indent=2))
    else:
        parser.error('Specify --refine or --verify')
