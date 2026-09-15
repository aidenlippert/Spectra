"""Build molecular Schur intervals using selected actions and occupation proof trees."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle, GapFailure
from experiments.marginal_molecular_gershgorin import schur_pivots, suggest_lower, SOURCES
from experiments.marginal_implicit_certificate import rational_text

ROOT = Path(__file__).resolve().parents[1]


def replay(certificate):
    if certificate.get('kind') != 'sparse_determinant_schur_v1':
        raise ValueError('Unsupported sparse determinant certificate')
    oracle = DeterminantOracle(certificate)
    p = certificate.get('retained_states')
    oracle.retained(p)
    gamma, lower = F(certificate['complement_lower']), F(certificate['lower'])
    if gamma <= lower:
        raise ValueError('Positive complementary Schur denominator required')
    tree = certificate.get('complement_tree')
    if type(tree) is not list:
        raise ValueError('Explicit complement proof tree required')
    _, stats = oracle.cover(p, gamma, tree)
    pivots = schur_pivots(oracle.retained_data(p, gamma), lower)
    if pivots is None:
        raise ValueError('Exact retained Schur positivity failed')
    upper = oracle.upper(certificate.get('independent_upper'))
    if lower > upper:
        raise ValueError('Inconsistent sparse determinant interval')
    return {'lower': rational_text(lower), 'upper': rational_text(upper),
            'width': rational_text(upper - lower), 'width_float': float(upper - lower),
            'retained_dimension': len(p), 'sector_dimension': math.comb(oracle.modes, oracle.particles),
            'complement_lower': rational_text(gamma), 'positive_denominator': rational_text(gamma - lower),
            'schur_pivots': [rational_text(x) for x in pivots], 'complement_coverage': stats,
            'unique_action_states': len(oracle.cache),
            'referenced_determinants': oracle.referenced_state_count(),
            'scope': 'Selected sparse determinant actions and complete occupation-tree complement proof; no full sector list or matrix; worst-case tree/action growth remains unbounded by a polynomial'}


def run_source(source_path, out, target=F(1, 10**7), selector='counterexample'):
    import numpy as np
    source_path, out = Path(source_path), Path(out)
    if selector not in ('counterexample', 'residual'):
        raise ValueError('Unknown sparse reference selector')
    if out.exists():
        raise ValueError('Preserve previous sparse molecular export')
    source = json.loads(source_path.read_text())
    # Only Hamiltonian and sector labels are imported; no source P or upper vector.
    base = {key: source[key] for key in ('modes', 'particles', 'hamiltonian')}
    oracle = DeterminantOracle(base)
    p, history, best = [(1 << oracle.particles) - 1], [], None
    started = time.monotonic()
    while True:
        block, leakage, _ = oracle.retained_data(p, F(0))
        _, vectors = np.linalg.eigh(np.array(block, dtype=float))
        amplitudes = [int(round(x * 10**10)) for x in vectors[:, 0]]
        witness = {'states': list(p), 'amplitudes': amplitudes}
        u = oracle.upper(witness)
        gamma = F(math.ceil(u * 10**10), 10**10) + F(1, 100)
        item = {'retained_dimension': len(p), 'upper': rational_text(u), 'gamma': rational_text(gamma)}
        try:
            tree, stats = oracle.cover(p, gamma)
        except GapFailure as failure:
            chosen = failure.state
            item.update(status='complement_counterexample', next_state=chosen, row_lower=rational_text(failure.lower))
            if selector == 'residual':
                residual = {}
                for state, amplitude in zip(p, amplitudes):
                    for destination, value in oracle.action(state).items():
                        if destination not in p:
                            residual[destination] = residual.get(destination, F(0)) + amplitude * value
                if residual:
                    candidate = max(residual, key=lambda s: (abs(residual[s]), -s))
                    # A proposal heuristic only; every final gap is still replayed.
                    if residual[candidate] ** 2 > F(1, 10**16) * sum(x * x for x in amplitudes):
                        chosen = candidate
                        item.update(status='ritz_residual', next_state=chosen, deferred_counterexample=failure.state)
        else:
            data = block, leakage, gamma
            lower = suggest_lower(data, u)
            item.update(status='accepted', lower=rational_text(lower), width_float=float(u - lower), coverage=stats)
            c = dict(base, kind='sparse_determinant_schur_v1', retained_states=list(p),
                     independent_upper=witness, complement_lower=rational_text(gamma),
                     lower=rational_text(lower), complement_tree=tree)
            if best is None or u - lower < F(best[1]['width']):
                best = c, replay(c)
            if u - lower <= target or len(p) == 32:
                history.append(item)
                break
            effective = np.array(block, dtype=float) - np.array(leakage, dtype=float) / float(gamma - lower)
            _, effective_vectors = np.linalg.eigh(effective)
            scores = {}
            for i, state in enumerate(p):
                for destination, value in oracle.action(state).items():
                    if destination not in p:
                        scores[destination] = scores.get(destination, 0.) + float(value) * effective_vectors[i, 0]
            if not scores:
                raise ValueError('No external coupling to enrich the retained reference')
            chosen = max(scores, key=lambda s: (abs(scores[s]), -s))
            item['next_state'] = chosen
        history.append(item)
        print(json.dumps({k: item[k] for k in ('retained_dimension', 'status')}), flush=True)
        if len(p) == 32:
            break
        p.append(chosen)
    if best is None:
        out.mkdir(parents=True)
        failure = {'status': 'not_accepted', 'reason': 'Retained budget exhausted without a complete gap certificate',
                   'retained_dimension': len(p), 'construction_unique_action_states': len(oracle.cache),
                   'construction_referenced_determinants': oracle.referenced_state_count(),
                   'retained_states': list(p), 'independent_upper': witness,
                   'last_gap_threshold': rational_text(gamma),
                   'elapsed_seconds': time.monotonic() - started, 'selector': selector,
                   'source_hamiltonian': str(source_path)}
        (out / 'failure.json').write_text(json.dumps(failure, indent=2) + '\n')
        (out / 'history.json').write_text(json.dumps(history, indent=2) + '\n')
        print(json.dumps(failure), flush=True)
        return
    certificate, receipt = best
    receipt.update(elapsed_seconds=time.monotonic() - started, construction_unique_action_states=len(oracle.cache),
                   construction_referenced_determinants=oracle.referenced_state_count(),
                   source_hamiltonian=str(source_path), target_width=rational_text(target), selector=selector,
                   stopping_reason='target_width' if F(receipt['width']) <= target else 'retained_dimension_budget')
    out.mkdir(parents=True)
    for filename, value in (('certificate.json', certificate), ('receipt.json', receipt), ('history.json', history)):
        (out / filename).write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ('width_float', 'retained_dimension', 'unique_action_states',
                                            'construction_unique_action_states', 'elapsed_seconds', 'stopping_reason')}), flush=True)


def run(name, target=F(1, 10**7), selector='counterexample'):
    if name not in SOURCES:
        raise ValueError('Unknown molecular fixture')
    out = ROOT / 'results/marginal_sparse_molecular' / (name + ('_residual' if selector == 'residual' else ''))
    return run_source(ROOT / SOURCES[name], out, target, selector)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=tuple(SOURCES), default='rectangle')
    parser.add_argument('--verify')
    parser.add_argument('--selector', choices=('counterexample', 'residual'), default='counterexample')
    parser.add_argument('--source')
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.source:
        if not args.output:
            parser.error('--source requires --output')
        run_source(args.source, args.output, selector=args.selector)
    else:
        if args.output:
            parser.error('--output requires --source')
        run(args.fixture, selector=args.selector)
