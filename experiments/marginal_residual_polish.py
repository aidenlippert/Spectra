"""Propose a scaled-DD residual and export it as exact rational pair atoms.

The original metric stays fixed. A new positive vector scales only the
residual comparison problem; acceptance uses the existing full exact replay.
"""
from fractions import Fraction as F
import json
from pathlib import Path

from experiments.marginal_signed_atoms import atom_vector, encode_vector, residual, verify_block, replay
from experiments.marginal_implicit_certificate import rational_text


def propose_residual_scaling(r, metric):
    """Numerical positive scaling proposal; the original metric is fixed."""
    import numpy as np
    n = len(r)
    rf = np.array(r, dtype=float)
    root_metric = np.sqrt(np.array(metric, dtype=float))
    comparison = -np.abs(rf)
    np.fill_diagonal(comparison, np.diag(rf))
    comparison /= root_metric[:, None] * root_metric[None, :]
    # Each connected residual component gets a strictly positive proposal.
    unseen = set(range(n)); scaling = [1] * n; numerical_floors = []
    while unseen:
        first = min(unseen); unseen.remove(first); group = [first]; todo = [first]
        while todo:
            i = todo.pop()
            adjacent = [j for j in sorted(unseen) if r[i][j]]
            unseen.difference_update(adjacent); group.extend(adjacent); todo.extend(adjacent)
        group.sort()
        eig, vectors = np.linalg.eigh(comparison[np.ix_(group, group)])
        numerical_floors.append(float(eig[0]))
        w = np.abs(vectors[:, 0])/root_metric[group]
        w /= np.max(w)
        for i, value in zip(group, w): scaling[i] = max(1, int(round(float(value)*10**9)))
    return scaling, min(numerical_floors)


def polish_block(a, item):
    r, metric = residual(a, item['metric_weights'], item['atoms'], item['atom_scale'], True)
    n = len(r)
    prior = min((r[i][i] - sum(abs(r[i][j]) for j in range(n) if j != i))/metric[i] for i in range(n))
    scaling, numerical_floor = propose_residual_scaling(r, metric)
    coefficients = {}
    for atom in item['atoms']:
        key = atom_vector(atom, n, True)
        coefficients[key] = coefficients.get(key, F(0)) + F(atom['weight'], item['atom_scale'])
    added = 0
    for i in range(n):
        for j in range(i+1, n):
            if not r[i][j]: continue
            wi, wj = scaling[i], scaling[j]
            maximum = max(wi, wj)
            sign = 1 if r[i][j] > 0 else -1
            key = ((i, j), (F(wj, maximum), F(sign*wi, maximum)))
            weight = abs(r[i][j])*F(maximum**2, wi*wj)
            coefficients[key] = coefficients.get(key, F(0)) + weight
            added += 1
    scale = 10**18
    atoms = []
    for key, value in coefficients.items():
        weight = (value*scale).numerator//(value*scale).denominator
        if weight: atoms.append(dict(encode_vector(key), weight=weight))
    candidate = dict(metric_weights=item['metric_weights'], atom_scale=scale, atoms=atoms)
    rr, mm = residual(a, candidate['metric_weights'], atoms, scale, True)
    proposed = min((rr[i][i] - sum(abs(rr[i][j]) for j in range(n) if j != i))/mm[i] for i in range(n))
    accepted = proposed >= prior
    result, bound = (candidate, proposed) if accepted else (item, prior)
    verify_block(a, result, bound, True)
    return result, bound, {'dimension': n, 'prior_exact_lower': rational_text(prior),
        'numerical_scaled_residual_lower': numerical_floor, 'scaling_weights': scaling,
        'pair_edges_processed': added, 'proposed_exact_lower': rational_text(proposed),
        'accepted': accepted, 'exported_exact_lower': rational_text(bound),
        'exported_exact_lower_float': float(bound), 'exported_atoms': len(result['atoms']),
        'scope': 'One scaled-DD residual proposal with the original metric fixed. Rational pair atoms are merged and quantized, then the complete exact residual is checked. No FW4 optimum or general scaling claim.'}


def construct(source, output, target=None):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    source, output = Path(source), Path(output)
    if output.exists(): raise ValueError('Preserve previous residual polish export')
    old = json.loads(source.read_text())
    replay(old)
    oracle = SpinZeroOracle(old)
    matrices = checked_blocks(oracle, old['retained_states'], old['blocks'])
    blocks = []; bounds = []; proposals = []
    for a, block in zip(matrices, old['blocks']):
        item, bound, diagnostic = polish_block(a, block)
        blocks.append(dict(item, states=block['states'])); bounds.append(bound); proposals.append(diagnostic)
    c = {k: old[k] for k in ('modes', 'particles', 'hamiltonian', 'retained_states')}
    c.update(kind='spin_rational_atom_complement_v1', blocks=blocks, target_lower=rational_text(min(bounds)))
    receipt = replay(c)
    receipt.update(source=str(source), source_target_lower=old['target_lower'])
    if target is not None:
        receipt.update(requested_target_lower=rational_text(F(target)), requested_target_certified=min(bounds)>=F(target))
    output.mkdir(parents=True)
    for name, data in [('certificate.json', c), ('receipt.json', receipt), ('proposal.json', proposals)]:
        (output/name).write_text(json.dumps(data, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True); parser.add_argument('--output', required=True)
    parser.add_argument('--target')
    args = parser.parse_args()
    print(json.dumps(construct(args.source, args.output, args.target), indent=2))
