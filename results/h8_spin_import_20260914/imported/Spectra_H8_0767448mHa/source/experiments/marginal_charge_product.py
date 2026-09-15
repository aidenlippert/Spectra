"""Positive charge-product metrics for complete weighted-DD complement proofs."""
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path


def feature_orbits(sites, family, max_pair_distance=None):
    if type(sites) is not int or not 2 <= sites <= 32 or family not in ('pair', 'pair_square'):
        raise ValueError('Supported bounded charge feature family required')
    distance = sites-1 if max_pair_distance is None else max_pair_distance
    if type(distance) is not int or not 1 <= distance < sites:
        raise ValueError('Pair distance must lie between one and sites minus one')
    labels = [((i, p),) for i in range(sites) for p in (1, 2)]
    labels += [((i, 1), (j, 1)) for i, j in combinations(range(sites), 2) if j-i <= distance]
    if family == 'pair_square':
        labels += [((i, 2), (j, 2)) for i, j in combinations(range(sites), 2) if j-i <= distance]
    groups = {}
    for label in labels:
        reflection = tuple(sorted((sites-1-i, p) for i, p in label))
        groups.setdefault(min(label, reflection), []).append(label)
    return list(groups.values())


def features(state, sites, family, max_pair_distance=None):
    q = [((state >> (2*i)) & 1) + ((state >> (2*i+1)) & 1) - 1 for i in range(sites)]
    result = []
    for orbit in feature_orbits(sites, family, max_pair_distance):
        value = 0
        for label in orbit:
            term = 1
            for i, power in label:
                term *= q[i]**power
            value += term
        result.append(value)
    return result


def expand_recipe(recipe, modes):
    if type(recipe) is not dict or recipe.get('kind') != 'spin_charge_product_complement_v1':
        raise ValueError('Expected charge-product complement recipe')
    if type(modes) is not int or modes % 4 or not 4 <= modes <= 64:
        raise ValueError('Even-site spin-orbital system required')
    rule = recipe['metric_rule']; family = rule['family']
    distance = rule.get('max_pair_distance')
    count = len(feature_orbits(modes//2, family, distance))
    raw = rule['factors']
    if type(raw) is not list or len(raw) != count or any(type(v) is not str for v in raw):
        raise ValueError('One exact factor per charge feature required')
    factors = [F(v) for v in raw]
    if any(v <= 0 for v in factors):
        raise ValueError('Strictly positive charge metric factors required')
    integer_scale = rule['integer_scale']
    if type(integer_scale) is not int or not 1 <= integer_scale <= 10**12:
        raise ValueError('Bounded positive integer metric scale required')
    if type(recipe['blocks']) is not list or not 1 <= len(recipe['blocks']) <= 384:
        raise ValueError('Bounded explicit complement blocks required')
    blocks = []
    for block in recipe['blocks']:
        if set(block) != {'states'} or type(block['states']) is not list or not 1 <= len(block['states']) <= 384:
            raise ValueError('Charge-product blocks contain states only')
        weights = []
        for state in block['states']:
            if type(state) is not int or not 0 <= state < 1 << modes:
                raise ValueError('Invalid determinant label')
            weight = F(1)
            for factor, power in zip(factors, features(state, modes//2, family, distance)):
                weight *= factor**power
            weights.append(weight)
        # Exact rational nearest-integer rounding; the existing metric gate
        # checks positivity and size, and replay verifies the rounded metric.
        blocks.append({'states': block['states'], 'metric_weights': [round(w*integer_scale) for w in weights],
                       'atoms': [], 'atom_scale': 1})
    return {'kind': 'spin_rational_atom_complement_v1', 'blocks': blocks}


def replay(certificate):
    from experiments.marginal_signed_atoms import replay as replay_atoms
    recipe = {k: certificate[k] for k in ('kind', 'metric_rule', 'blocks')}
    expanded = expand_recipe(recipe, certificate['modes'])
    result = replay_atoms(dict(certificate, **expanded))
    result.update(metric_factor_count=len(recipe['metric_rule']['factors']),
                  metric_family=recipe['metric_rule']['family'],
                  metric_scope='Exactly rounded product of reflection-tied positive one-site and two-site charge factors. Reflection ties metric parameters only; H symmetry is not assumed. All complement rows are still explicitly rebuilt and the actual integer metric is verified.')
    return result


def propose(source, output, family='pair_square', target=F('-6.24'), denominator=10**8, *, max_pair_distance=None):
    import numpy as np
    from scipy.optimize import minimize
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError('Preserve previous charge metric export')
    if type(denominator) is not int or denominator < 1:
        raise ValueError('Positive factor denominator required')
    c = json.loads(source.read_text()); oracle = SpinZeroOracle(c)
    if len(c['blocks']) != 1:
        raise ValueError('Proposer currently uses one complete Q block')
    matrices = checked_blocks(oracle, c['retained_states'], c['blocks'], max_block_dimension=384)
    states = c['blocks'][0]['states']
    x = np.array([[float(v) for v in row] for row in matrices[0]])
    diagonal = x.diagonal().copy(); off = np.abs(x-np.diag(diagonal))
    f = np.array([features(s, oracle.modes//2, family, max_pair_distance) for s in states], dtype=float)
    # Zero initialization uses no inherited metric or spectral vector.
    def margins(z):
        w = np.exp(f @ z[:-1])
        return diagonal-(off @ w)/w-z[-1]
    def jacobian(z):
        w = np.exp(f @ z[:-1]); radius = (off @ w)/w
        return np.column_stack([radius[:, None]*f-(off @ (w[:, None]*f))/w[:, None], -np.ones(len(states))])
    start = np.r_[np.zeros(f.shape[1]), min(diagonal-off.sum(axis=1))]
    result = minimize(lambda z: -z[-1], start, jac=lambda z: np.r_[np.zeros(f.shape[1]), -1],
                      constraints=[{'type': 'ineq', 'fun': margins, 'jac': jacobian}],
                      bounds=[(-5, 5)]*f.shape[1]+[(None, None)], method='SLSQP',
                      options={'maxiter': 400, 'ftol': 1e-11})
    if not result.success:
        raise ValueError('No successful metric proposal: '+result.message)
    factors = [F(round(float(np.exp(v))*denominator), denominator) for v in result.x[:-1]]
    cert = {k: c[k] for k in ('modes', 'particles', 'hamiltonian', 'retained_states')}
    cert.update(kind='spin_charge_product_complement_v1', target_lower=str(F(target)),
                blocks=[{'states': states}], metric_rule={'family': family, 'factors': [str(v) for v in factors], 'integer_scale': 10**8})
    if max_pair_distance is not None:
        cert['metric_rule']['max_pair_distance'] = max_pair_distance
    receipt = replay(cert)
    output.mkdir(parents=True)
    for name, data in [('certificate', cert), ('receipt', receipt),
                       ('proposal', {'numerical_lower': float(min(margins(np.r_[result.x[:-1], 0]))),
                                     'optimizer_success': bool(result.success), 'feature_count': f.shape[1],
                                     'initialization': 'Zero feature parameters; no inherited metric.',
                                     'scope': 'Numerical proposal, exact replay determines acceptance; no optimality claim.'})]:
        (output/(name+'.json')).write_text(json.dumps(data, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('--verify'); p.add_argument('--source'); p.add_argument('--output')
    p.add_argument('--family', default='pair_square'); p.add_argument('--target', default='-6.24')
    p.add_argument('--max-pair-distance', type=int)
    a = p.parse_args()
    print(json.dumps(replay(json.loads(Path(a.verify).read_text())) if a.verify
                     else propose(a.source, a.output, a.family, F(a.target), max_pair_distance=a.max_pair_distance), indent=2))
