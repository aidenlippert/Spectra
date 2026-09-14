"""Construct a complete spin-reduced interval directly from a bare Hamiltonian."""
from fractions import Fraction as F
from itertools import combinations
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_h6_complement import checked_blocks, factor_blocks
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_sparse_response import q_action, prepare, response_matrix, refine_response
from experiments.marginal_sparse_upper import refine as refine_upper
from experiments.marginal_spin_reduction import SpinZeroOracle, symmetrize, replay
from experiments.marginal_symbolic import encode, add, scale


def select_reference(oracle, cap, *, seed=None):
    import numpy as np
    if type(cap) is not int or not 1 <= cap <= 32:
        raise ValueError('Reference budget must be one to 32')
    initial = (1 << oracle.particles) - 1 if seed is None else seed
    oracle.retained([initial])
    p, history = [initial], []
    for _ in range(min(cap, oracle.sector_dimension - 1)):
        a = [[oracle.action(t).get(s, F(0)) for t in p] for s in p]
        _, vectors = np.linalg.eigh(np.array(a, dtype=float))
        integers = [int(round(x * 10**10)) for x in vectors[:, 0]]
        witness = {'states': list(p), 'amplitudes': integers}
        upper = oracle.upper(witness)
        residual = q_action(oracle, set(p), dict(zip(p, map(F, integers))))
        history.append({'retained_dimension': len(p), 'upper': rational_text(upper), 'external_support': len(residual)})
        if len(p) == cap or len(p) == oracle.sector_dimension - 1 or not residual:
            break
        p.append(min(residual, key=lambda s: (-abs(residual[s]), s)))
    return p, witness, history


def spin_states(oracle):
    if oracle.sector_dimension > 4096:
        raise ValueError('Explicit spin-sector construction exceeds its 4096-state budget')
    occupied = [sum(1 << (2*i) for i in group)
                for group in combinations(range(oracle.modes // 2), oracle.particles // 2)]
    return sorted(a | (b << 1) for a in occupied for b in occupied)


def spin_blocks(oracle, p):
    remaining = set(spin_states(oracle)) - oracle.retained(p)
    groups = []
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        queue, group = [seed], []
        while queue:
            state = queue.pop()
            group.append(state)
            adjacent = sorted(s for s in oracle.action(state) if s in remaining)
            remaining.difference_update(adjacent)
            queue.extend(adjacent)
        groups.append(sorted(group))
    return groups


def wrap(original, reduced):
    return dict(original, kind='spin_reduced_perturbation_interval_v1', spin_symmetric_certificate=reduced)


def build(source, output, *, reference_cap=32, upper_steps=20, target=F(1, 10**7)):
    source, out = Path(source), Path(output)
    if out.exists():
        raise ValueError('Preserve previous direct spin construction')
    # Deliberately import only these three fields, never a supplied P, upper, gap, or factor.
    input_data = json.loads(source.read_text())
    original = {key: input_data[key] for key in ('modes', 'particles', 'hamiltonian')}
    original_oracle = DeterminantOracle(original)
    symmetric = dict(original, hamiltonian=encode(symmetrize(original_oracle.h, original_oracle.modes)))
    oracle = SpinZeroOracle(symmetric)
    if oracle.sector_dimension > 4096:
        raise ValueError('Explicit spin-sector construction exceeds its 4096-state budget')
    error = sum(abs(x) for x in add(original_oracle.h, scale(oracle.h, -1)).values())
    if target <= 2 * error:
        raise ValueError('Requested width is below the conservative spin-transfer budget')
    out.mkdir(parents=True)
    started = time.monotonic()
    try:
        _build(source, out, original, symmetric, oracle, error, reference_cap, upper_steps, F(target), started)
    except Exception as exc:
        (out / 'failure.json').write_text(json.dumps({'status': 'not_accepted', 'error': str(exc),
            'elapsed_seconds': time.monotonic() - started, 'input': str(source)}, indent=2) + '\n')
        raise


def _build(source, out, original, symmetric, oracle, error, cap, steps, target, started):
    import numpy as np
    p, witness, selection = select_reference(oracle, cap)
    (out / 'symmetric_hamiltonian.json').write_text(json.dumps(symmetric, indent=2) + '\n')
    (out / 'selection.json').write_text(json.dumps({'retained_states': p, 'independent_upper': witness,
        'history': selection, 'unique_action_states': len(oracle.cache),
        'referenced_determinants': oracle.referenced_state_count()}, indent=2) + '\n')
    refine_upper(out / 'symmetric_hamiltonian.json', out / 'selection.json', out / 'upper', steps)
    upper_certificate = json.loads((out / 'upper/certificate.json').read_text())
    upper = oracle.upper(upper_certificate['independent_upper'])
    groups = spin_blocks(oracle, p)
    matrices = checked_blocks(oracle, p, [{'states': group} for group in groups])
    minimum = min(float(np.linalg.eigvalsh(np.array(a, dtype=float))[0]) for a in matrices)
    gamma = F(math.floor((minimum - 1e-4) * 10**12), 10**12)
    if gamma <= upper:
        raise ValueError('Selected reference has no certified complement threshold above the upper')
    blocks = factor_blocks(groups, matrices, gamma)
    c = dict(symmetric, kind='spin_zero_response_recipe_v1', retained_states=p,
             independent_upper=upper_certificate['independent_upper'], complement_lower=rational_text(gamma), blocks=blocks)
    data = prepare(oracle, p, gamma, [])
    for exponent in range(32):
        lower = F(math.floor(min(upper, gamma) - 2**exponent))
        if ldl_pivots(response_matrix(data, lower)) is not None:
            break
    else:
        raise ValueError('Initial exact scalar Schur bound failed')
    c['lower'] = rational_text(lower)
    def check(candidate, *, audit=None):
        return replay(wrap(symmetric, candidate), audit=audit)
    refine_response(c, out / 'response_recipe', {'upper': rational_text(upper), 'width_float': float(upper - lower)},
                    target - 2*error, str(source), oracle=oracle, check=check)
    reduced = json.loads((out / 'response_recipe/certificate.json').read_text())
    final = wrap(original, reduced)
    receipt = replay(final)
    if F(receipt['width']) > target:
        raise ValueError('Direct construction exhausted its response budget before the original-H target')
    all_states = set(spin_states(oracle))
    if len(oracle.cache) != oracle.sector_dimension or not set(oracle.cache) <= all_states:
        raise ValueError('Unexpected direct-construction state coverage')
    receipt.update(input=str(source), imported_fields=['modes', 'particles', 'hamiltonian'],
        elapsed_seconds=time.monotonic() - started, selection_dimension=len(p),
        upper_steps=steps, gamma_proposal_control=minimum,
        largest_numerical_eigensolve=max(len(p), max(map(len, groups)), steps),
        construction_unique_determinants=receipt['referenced_determinants'],
        construction_unique_source_determinants=len(all_states),
        construction_spin_symmetric_sources=oracle.sector_dimension,
        construction_original_sources=receipt['original_action_states'],
        timing_scope='Timer starts after input validation and spin projection; repeated exact replays are included.',
        discovery_scope='Complete construction from bare H and sector labels: spin projection, residual-selected P, own Krylov upper, spin-zero Q graph, newly proposed factors and response. Enumerates all spin-zero determinants, never a full fixed-N list; original-H witness actions may reference additional states when Sz is broken. No inherited proof, FCI input, or asymptotic efficiency claim.')
    (out / 'certificate.json').write_text(json.dumps(final, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ('width_float', 'response_dimension', 'construction_unique_determinants', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference-cap', type=int, default=32)
    parser.add_argument('--upper-steps', type=int, default=20)
    args = parser.parse_args()
    build(args.source, args.output, reference_cap=args.reference_cap, upper_steps=args.upper_steps)
