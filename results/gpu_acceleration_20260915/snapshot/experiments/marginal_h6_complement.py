"""Comparison-family obstruction and sign-preserving finite H6 complement proof."""
from fractions import Fraction as F
from itertools import combinations
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_fixed_point_ldl import propose, verify
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_sparse_response import parse_basis

ROOT = Path(__file__).resolve().parents[1]


def comparison_rayleigh(oracle, retained, vector):
    energy = F(0)
    for s, x in vector.items():
        for t, value in oracle.action(s).items():
            if t not in retained and t in vector:
                energy += x * vector[t] * (value if t == s else -abs(value))
    return energy / sum(x*x for x in vector.values())


def replay_comparison(certificate):
    if certificate.get('kind') != 'comparison_gershgorin_obstruction_v1':
        raise ValueError('Unsupported comparison obstruction')
    oracle = DeterminantOracle(certificate)
    retained = oracle.retained(certificate.get('retained_states'))
    vector = parse_basis(oracle, retained, [certificate.get('comparison_witness')])[0]
    target = F(certificate['target_lower'])
    value = comparison_rayleigh(oracle, retained, vector)
    if value >= target:
        raise ValueError('Comparison witness does not obstruct the target')
    return {'target_lower': rational_text(target), 'comparison_rayleigh': rational_text(value),
            'deficit': rational_text(target - value), 'deficit_float': float(target - value),
            'witness_support': len(vector), 'unique_action_states': len(oracle.cache),
            'referenced_determinants': oracle.referenced_state_count(),
            'scope': 'No positive diagonal weighting of determinant-basis Gershgorin can reach this target for the specified Q; no obstruction to the physical spectral gap is claimed'}


def checked_blocks(oracle, p, blocks, *, max_block_dimension=256):
    if type(max_block_dimension) is not int or not 1 <= max_block_dimension <= 384:
        raise ValueError('Bounded complement dimension required')
    retained = oracle.retained(p)
    if type(blocks) is not list or not blocks:
        raise ValueError('Complete complement block recipe required')
    groups = []
    for item in blocks:
        if type(item) is not dict or type(item.get('states')) is not list or not 1 <= len(item['states']) <= max_block_dimension:
            raise ValueError('Invalid bounded complement block')
        states = item['states']
        if any(not oracle.valid_state(s) or s in retained for s in states):
            raise ValueError('Complement blocks must contain physical Q determinants')
        groups.append(states)
    flattened = [s for group in groups for s in group]
    if len(flattened) != oracle.sector_dimension - len(retained) or len(set(flattened)) != len(flattened):
        raise ValueError('Blocks must cover every complementary determinant exactly once')
    membership = {s: i for i, group in enumerate(groups) for s in group}
    matrices = []
    for i, group in enumerate(groups):
        columns = [oracle.action(s) for s in group]
        if any(t not in retained and membership[t] != i for column in columns for t in column):
            raise ValueError('Declared complement blocks have a nonzero interblock coupling')
        matrices.append([[columns[j].get(s, F(0)) for j in range(len(group))] for s in group])
    return matrices


def replay_factor(certificate):
    if certificate.get('kind') != 'factor_complement_bound_v1':
        raise ValueError('Unsupported complement factor proof')
    oracle = DeterminantOracle(certificate)
    gamma = F(certificate['complement_lower'])
    result = verify_complement(oracle, certificate.get('retained_states'), gamma, certificate.get('blocks'))
    result.update(unique_action_states=len(oracle.cache), referenced_determinants=oracle.referenced_state_count())
    return result


def verify_complement(oracle, p, gamma, blocks):
    """Bind the complete factor proof directly to the caller's Hamiltonian and P."""
    matrices = checked_blocks(oracle, p, blocks)
    receipts = []
    for a, item in zip(matrices, blocks):
        shifted = [[x - gamma * (i == j) for j, x in enumerate(row)] for i, row in enumerate(a)]
        r = verify(shifted, item.get('factor'))
        receipts.append({key: rational_text(value) if isinstance(value, F) else value for key, value in r.items()})
    return {'complement_lower': rational_text(gamma), 'block_dimensions': [len(a) for a in matrices],
            'block_receipts': receipts,
            'scope': 'Exact positive complement gap from complete coordinate-block coverage and verified fixed-point LDL residuals; all Q states are explicitly represented; not a ground-energy interval or scaling result'}


def factor_blocks(groups, matrices, gamma):
    blocks = []
    for group, a in zip(groups, matrices):
        shifted = [[x - gamma * (i == j) for j, x in enumerate(row)] for i, row in enumerate(a)]
        denominator = math.lcm(*(F(x).denominator for row in shifted for x in row))
        for digits in (12, 24, 40, 56):
            try:
                factor = propose(shifted, math.lcm(10**digits, denominator))
                margin = verify(shifted, factor)
            except ValueError:
                continue
            blocks.append({'states': group, 'factor': factor})
            print(json.dumps({'block_dimension': len(group), 'digits': digits, 'margin': float(margin['margin'])}), flush=True)
            break
        else:
            raise ValueError('Fixed-point complement factor budget exhausted')
    return blocks


def retarget(path, output, gamma):
    source, out, gamma = Path(path), Path(output), F(gamma)
    if out.exists():
        raise ValueError('Preserve previous complement export')
    c = json.loads(source.read_text())
    if c.get('kind') != 'factor_complement_bound_v1':
        raise ValueError('Expected a complement factor source')
    oracle = DeterminantOracle(c)
    started = time.monotonic()
    matrices = checked_blocks(oracle, c.get('retained_states'), c.get('blocks'))
    groups = [item['states'] for item in c['blocks']]
    c.update(complement_lower=rational_text(gamma), blocks=factor_blocks(groups, matrices, gamma))
    r = replay_factor(c)
    r.update(source=str(source), elapsed_seconds=time.monotonic() - started)
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(c, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(r, indent=2) + '\n')
    print(json.dumps({'gamma': float(gamma), 'elapsed_seconds': r['elapsed_seconds']}), flush=True)


def run():
    import numpy as np
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import eigsh
    root = ROOT / 'results/marginal_h6'
    out = root / 'complement_proofs'
    if out.exists():
        raise ValueError('Preserve previous complement proofs')
    fixture = json.loads((root / 'fixture.json').read_text())
    p = json.loads((root / 'sparse_p32/failure.json').read_text())['retained_states']
    upper = F(json.loads((root / 'krylov_upper_20/receipt.json').read_text())['upper'])
    gamma = F(math.floor((upper + F(1, 100)) * 10**12), 10**12)
    base = {key: fixture[key] for key in ('hamiltonian', 'modes', 'particles')}
    base['retained_states'] = p
    oracle = DeterminantOracle(base)
    retained = oracle.retained(p)
    # Explicit enumeration is a finite diagnostic/control here, not the sparse constructor.
    q = sorted(sum(1 << i for i in group) for group in combinations(range(oracle.modes), oracle.particles))
    q = [s for s in q if s not in retained]
    index = {s: i for i, s in enumerate(q)}
    rows, columns, physical, comparison = [], [], [], []
    started = time.monotonic()
    for s, j in index.items():
        for t, value in oracle.action(s).items():
            if t in index:
                rows.append(index[t]); columns.append(j)
                physical.append(float(value)); comparison.append(float(value) if s == t else -abs(float(value)))
    a = coo_matrix((physical, (rows, columns)), shape=(len(q), len(q))).tocsr()
    c = coo_matrix((comparison, (rows, columns)), shape=a.shape).tocsr()
    ep, vp = eigsh(a, k=1, which='SA', tol=1e-12, v0=np.ones(len(q)))
    ec, vc = eigsh(c, k=1, which='SA', tol=1e-12, v0=np.ones(len(q)))
    ranked = sorted(range(len(q)), key=lambda i: (-abs(vc[i, 0]), q[i]))
    obstruction = None
    for count in (4, 8, 12, 16, 24, 32, 64, 128, 256, len(q)):
        chosen = sorted(ranked[:count], key=lambda i: q[i])
        witness = {'states': [q[i] for i in chosen],
                   'amplitudes': [int(round(abs(vc[i, 0]) * 10**10)) for i in chosen]}
        trial = dict(base, kind='comparison_gershgorin_obstruction_v1', comparison_witness=witness, target_lower=rational_text(gamma))
        try:
            receipt = replay_comparison(trial)
        except ValueError:
            continue
        obstruction = trial, receipt
        break
    if obstruction is None:
        raise ValueError('No comparison obstruction exported')
    remaining, groups = set(q), []
    while remaining:
        seed = min(remaining)
        group, queue = [], [seed]
        remaining.remove(seed)
        while queue:
            state = queue.pop()
            group.append(state)
            adjacent = sorted(t for t in oracle.action(state) if t in remaining)
            remaining.difference_update(adjacent)
            queue.extend(adjacent)
        groups.append(sorted(group))
    matrices = checked_blocks(oracle, p, [{'states': group} for group in groups])
    blocks = factor_blocks(groups, matrices, gamma)
    certificate = dict(base, kind='factor_complement_bound_v1', complement_lower=rational_text(gamma), blocks=blocks)
    receipt = replay_factor(certificate)
    diagnostic = {'status': 'numerical_full_Q_control_with_separate_exact_proofs', 'q_dimension': len(q),
                  'physical_q_minimum_numerical': float(ep[0]), 'comparison_q_minimum_numerical': float(ec[0]),
                  'physical_residual_norm': float(np.linalg.norm(a @ vp[:, 0] - ep[0] * vp[:, 0])),
                  'comparison_residual_norm': float(np.linalg.norm(c @ vc[:, 0] - ec[0] * vc[:, 0])),
                  'elapsed_seconds': time.monotonic() - started,
                  'scope': 'Explicit 892-state sparse Q matrices used for numerical controls and obstruction discovery; numerical eigenvalues are not certified endpoints'}
    out.mkdir(parents=True)
    for filename, value in (('comparison_certificate.json', obstruction[0]), ('comparison_receipt.json', obstruction[1]),
                            ('factor_certificate.json', certificate), ('factor_receipt.json', receipt), ('diagnostic.json', diagnostic)):
        (out / filename).write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps({'gamma': float(gamma), 'comparison_support': obstruction[1]['witness_support'],
                      'comparison_deficit': obstruction[1]['deficit_float'], 'elapsed_seconds': diagnostic['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify-comparison')
    parser.add_argument('--verify-factor')
    parser.add_argument('--retarget')
    parser.add_argument('--output')
    parser.add_argument('--gamma')
    args = parser.parse_args()
    if args.verify_comparison:
        print(json.dumps(replay_comparison(json.loads(Path(args.verify_comparison).read_text())), indent=2))
    elif args.verify_factor:
        print(json.dumps(replay_factor(json.loads(Path(args.verify_factor).read_text())), indent=2))
    elif args.retarget and args.output and args.gamma:
        retarget(args.retarget, args.output, args.gamma)
    else:
        run()
