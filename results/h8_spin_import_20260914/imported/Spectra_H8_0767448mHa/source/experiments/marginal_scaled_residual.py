"""Compact exact positivity from sparse atoms and a positive residual scaling.

Each residual edge's rational pair factor is implicit in one positive vector.
Replay uses only standard-library arithmetic and full physical Q coverage.
"""
from fractions import Fraction as F
import json
from pathlib import Path

from experiments.marginal_signed_atoms import residual, replay as replay_atoms
from experiments.marginal_implicit_certificate import rational_text

MAX_SCALING = 10**12


def scaled_bound(r, metric, scaling):
    n = len(r)
    if (type(scaling) is not list or len(scaling) != n
            or any(type(x) is not int or not 1 <= x <= MAX_SCALING for x in scaling)):
        raise ValueError('One bounded positive integer residual scaling per coordinate required')
    margins = [r[i][i] - sum(abs(r[i][j])*F(scaling[j], scaling[i])
                            for j in range(n) if j != i) for i in range(n)]
    return min(value/m for value, m in zip(margins, metric))


def refine_scaling(r, metric, scaling, steps=16):
    """Bounded exact improvements of the weakest scaled residual row.

    Common multiplication changes no ratio. Incrementing one coordinate
    after that multiplication allows finer rational steps without floats.
    Only strict improvements of the full minimum are retained.
    """
    if type(steps) is not int or not 0 <= steps <= 64:
        raise ValueError('Residual refinement budget must be zero through 64')
    bound = scaled_bound(r, metric, scaling)
    t = list(scaling); n = len(r); history = []
    cost = [sum((abs(r[i][j])*t[j] for j in range(n) if j != i), F(0)) for i in range(n)]
    for _ in range(steps):
        rows = [(r[i][i]-cost[i]/t[i])/metric[i] for i in range(n)]
        k = min(range(n), key=rows.__getitem__)
        best = None
        for factor in (1, 2, 10):
            candidate = [factor*x+int(i==k) for i,x in enumerate(t)]
            if max(candidate) > MAX_SCALING: continue
            updated = [factor*cost[i]+(abs(r[i][k]) if i != k else 0) for i in range(n)]
            value = min((r[i][i]-updated[i]/candidate[i])/metric[i] for i in range(n))
            if value > bound and (best is None or value > best[0]):
                best = value, candidate, updated, factor
        if best is None: break
        bound, t, cost, factor = best
        history.append({'row': k, 'common_factor': factor, 'exact_lower': rational_text(bound)})
    return t, bound, history


def verify_block(a, item, gamma):
    r, metric = residual(a, item.get('metric_weights'), item.get('atoms'), item.get('atom_scale'), True)
    bound = scaled_bound(r, metric, item.get('residual_scaling'))
    if bound < F(gamma): raise ValueError('Scaled residual does not certify the requested threshold')
    return {'dimension': len(r), 'explicit_atom_count': len(item['atoms']),
            'explicit_maximum_support': max([0]+[len(v['indices']) for v in item['atoms']]),
            'implicit_pair_count': sum(bool(r[i][j]) for i in range(len(r)) for j in range(i+1, len(r))),
            'scaled_residual_lower': rational_text(bound), 'scaled_residual_lower_float': float(bound),
            'minimum_margin': rational_text(bound-F(gamma))}


def replay(certificate):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    if certificate.get('kind') != 'spin_scaled_residual_complement_v1':
        raise ValueError('Unsupported scaled residual certificate')
    oracle = SpinZeroOracle(certificate)
    blocks = certificate.get('blocks')
    matrices = checked_blocks(oracle, certificate.get('retained_states'), blocks,max_block_dimension=384)
    gamma = F(certificate['target_lower'])
    rows = [verify_block(a, block, gamma) for a, block in zip(matrices, blocks)]
    return {'target_lower': rational_text(gamma), 'target_lower_float': float(gamma), 'blocks': rows,
            'unique_action_states': len(oracle.cache), 'referenced_determinants': oracle.referenced_state_count(),
            'target_certified': True,
            'scope': 'Exact complete-Q positivity from explicit rational atoms on at most four coordinates and a positive rational scaling of the full residual. Each residual edge defines an implicit positive two-coordinate factor. The original W² metric is fixed; no solver, spectral approximation, FW4 optimality, or general scaling claim is used in replay.'}


def construct(source, output, scaling_source=None, target=None, refine_steps=16):
    from experiments.marginal_spin_reduction import SpinZeroOracle
    from experiments.marginal_h6_complement import checked_blocks
    source, output = Path(source), Path(output)
    if output.exists(): raise ValueError('Preserve previous scaled residual export')
    old = json.loads(source.read_text())
    compact = old.get('kind') == 'spin_scaled_residual_complement_v1'
    (replay if compact else replay_atoms)(old)
    proposals = None
    if scaling_source is not None:
        proposals = json.loads(Path(scaling_source).read_text())
        if (type(proposals) is not list or len(proposals) != len(old['blocks'])
                or any(type(p) is not dict or p.get('dimension') != len(b['states'])
                       for p, b in zip(proposals, old['blocks']))):
            raise ValueError('One dimension-matched scaling proposal per source block required')
    oracle = SpinZeroOracle(old)
    matrices = checked_blocks(oracle, old['retained_states'], old['blocks'],max_block_dimension=384)
    blocks = []; bounds = []; history = []
    for index, (a, block) in enumerate(zip(matrices, old['blocks'])):
        r, metric = residual(a, block['metric_weights'], block['atoms'], block['atom_scale'], True)
        initial_scaling = block['residual_scaling'] if compact else [1]*len(a)
        initial = scaled_bound(r, metric, initial_scaling)
        numerical = None
        if proposals is None:
            from experiments.marginal_residual_polish import propose_residual_scaling
            scaling, numerical = propose_residual_scaling(r, metric)
        else:
            scaling = proposals[index].get('scaling_weights')
        scaling, bound, refinement = refine_scaling(r, metric, scaling, refine_steps)
        accepted = bound >= initial
        if not accepted: scaling, bound = initial_scaling, initial
        item = {key: block[key] for key in ('states', 'metric_weights', 'atom_scale', 'atoms')}
        item['residual_scaling'] = scaling
        verify_block(a, item, bound)
        blocks.append(item); bounds.append(bound)
        history.append({'dimension': len(a), 'initial_exact_lower': rational_text(initial),
                        'numerical_proposal_lower': numerical, 'accepted_scaling': accepted,
                        'exact_refinement': refinement,
                        'exported_exact_lower': rational_text(bound), 'exported_exact_lower_float': float(bound)})
    certificate = {key: old[key] for key in ('modes', 'particles', 'hamiltonian', 'retained_states')}
    certificate.update(kind='spin_scaled_residual_complement_v1', blocks=blocks, target_lower=rational_text(min(bounds)))
    receipt = replay(certificate)
    receipt.update(source=str(source), source_target_lower=old['target_lower'],
                   scaling_proposal_source=str(scaling_source) if scaling_source is not None else None)
    if target is not None:
        receipt.update(requested_target_lower=rational_text(F(target)), requested_target_certified=min(bounds)>=F(target))
    output.mkdir(parents=True)
    for name, data in [('certificate.json', certificate), ('receipt.json', receipt), ('proposal.json', history)]:
        (output/name).write_text(json.dumps(data, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify'); parser.add_argument('--source'); parser.add_argument('--output')
    parser.add_argument('--scaling-source'); parser.add_argument('--target'); parser.add_argument('--refine-steps', type=int, default=16)
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.source and args.output:
        print(json.dumps(construct(args.source, args.output, args.scaling_source, args.target, args.refine_steps), indent=2))
    else: parser.error('Specify --verify or --source and --output')
