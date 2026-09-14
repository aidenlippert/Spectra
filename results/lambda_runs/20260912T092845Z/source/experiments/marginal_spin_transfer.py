"""Connected spin-preserving transfer without factoring the joined Q block."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_sparse_response import prepare, response_matrix, refine_response
from experiments.marginal_sparse_upper import refine as refine_upper
from experiments.marginal_spin_constructor import wrap, spin_blocks
from experiments.marginal_spin_reduction import SpinZeroOracle, replay, spin_gap
from experiments.marginal_symbolic import add, mono, encode, scale
from experiments.marginal_implicit_certificate import rational_text


def transfer(source, output, strength=F(1, 1000), upper_steps=24, target=F(1, 10**7), *, norm_squares=False, temple_offset=None):
    source, out, strength = Path(source), Path(output), F(strength)
    if out.exists():
        raise ValueError('Preserve previous connected spin transfer')
    previous = json.loads(source.read_text())
    if previous.get('kind')=='spin_temple_interval_v1':
        from experiments.marginal_spin_temple import replay as temple_replay
        source_check=temple_replay(previous)
    else:
        source_check = replay(previous)
    base = previous['spin_symmetric_certificate']
    if 'complement_reference' in base or not any(key in base for key in ('blocks','complement_atoms')) or base['modes'] < 4:
        raise ValueError('A direct reference block certificate on at least four modes is required')
    delta = add(*(mono(((1, i), (0, j)), strength) for i, j in ((0, 2), (2, 0), (1, 3), (3, 1))))
    original_oracle, reference = DeterminantOracle(previous), SpinZeroOracle(base)
    original = {key: previous[key] for key in ('modes', 'particles')}
    original['hamiltonian'] = encode(add(original_oracle.h, delta))
    symmetric = dict(original, hamiltonian=encode(add(reference.h, delta)))
    oracle = SpinZeroOracle(symmetric)
    norm_certificate = None
    if norm_squares:
        from experiments.marginal_operator_norm import hopping_norm, verify_norm
        norm_certificate = hopping_norm(strength)
        eta, _ = verify_norm(delta, base['modes'], norm_certificate)
    else:
        eta = sum(abs(x) for x in delta.values())
    gamma = F(base['complement_lower']) - eta
    error = sum(abs(x) for x in add(original_oracle.h, scale(reference.h, -1)).values())
    if target <= 2*error:
        raise ValueError('Insufficient target budget for original-H transfer')
    out.mkdir(parents=True)
    started = time.monotonic()
    try:
        _transfer(source, out, base, original, symmetric, oracle, gamma, error, delta, upper_steps, F(target), started, source_check, norm_certificate, temple_offset)
    except Exception as exc:
        (out / 'failure.json').write_text(json.dumps({'status': 'not_accepted', 'error': str(exc)}, indent=2) + '\n')
        raise


def _transfer(source, out, base, original, symmetric, oracle, gamma, error, delta, steps, target, started, source_check, norm_certificate, temple_offset):
    (out / 'symmetric_hamiltonian.json').write_text(json.dumps(symmetric, indent=2) + '\n')
    (out / 'seed.json').write_text(json.dumps({'independent_upper': base['independent_upper']}, indent=2) + '\n')
    refine_upper(out / 'symmetric_hamiltonian.json', out / 'seed.json', out / 'upper', steps)
    upper_recipe = json.loads((out / 'upper/certificate.json').read_text())['independent_upper']
    upper = oracle.upper(upper_recipe)
    if gamma <= upper:
        raise ValueError('Transferred reference gap does not exceed the refined upper')
    p = base['retained_states']
    proof_key='complement_atoms' if 'complement_atoms' in base else 'blocks'
    c = dict(symmetric, kind='spin_zero_response_recipe_v1', retained_states=p, independent_upper=upper_recipe,
        complement_lower=rational_text(gamma), complement_reference={key: base[key] for key in ('hamiltonian', 'complement_lower', proof_key)})
    if norm_certificate is not None:
        c['complement_reference']['norm_certificate'] = norm_certificate
    _, gap_oracle = spin_gap(oracle, p, c)
    if temple_offset is not None:
        from experiments.marginal_spin_temple import discover
        final, r, history = discover(original, c, oracle, F(temple_offset), target)
        (out/'response_history.json').write_text(json.dumps(history, indent=2)+'\n')
    else:
        data = prepare(oracle, p, gamma, [])
        for exponent in range(32):
            lower = F(math.floor(min(upper, gamma) - 2**exponent))
            if ldl_pivots(response_matrix(data, lower)) is not None:
                break
        else:
            raise ValueError('Initial exact transferred Schur bound failed')
        c['lower'] = rational_text(lower)
        def check(candidate, *, audit=None):
            return replay(wrap(symmetric, candidate), audit=audit)
        refine_response(c, out / 'response_recipe', {'upper': rational_text(upper), 'width_float': float(upper-lower)},
                        target-2*error, str(source), oracle=oracle, check=check)
        reduced = json.loads((out / 'response_recipe/certificate.json').read_text())
        final = wrap(original, reduced)
        r = replay(final)
    if F(r['width']) > target:
        raise ValueError('Transferred interval exhausted the response budget before target')
    groups = spin_blocks(oracle, p)
    r.update(source=str(source), perturbation=encode(delta), elapsed_seconds=time.monotonic()-started,
        current_q_component_dimensions=list(map(len, groups)),
        largest_factor_dimension=max(len(item['states']) for item in base['blocks']) if 'blocks' in base else 0,
        construction_unique_source_determinants=len(set(oracle.cache) | set(gap_oracle.cache)),
        source_verification_action_states={'source_spin_symmetric': source_check['spin_symmetric_action_states'],
                                          'source_original': source_check['original_action_states']},
        construction_gap_reference_action_states=len(gap_oracle.cache),
        construction_current_spin_action_states=len(oracle.cache),
        construction_current_original_action_states=r['original_action_states'],
        discovery_scope='Transfer from a certified reference: old P and reference factors are reused; upper and response are regenerated on the connected Hamiltonian. No joined-Q factorization or full-sector eigensolve. Explicit spin-zero action coverage remains.',
        timing_scope='After source-certificate replay and input preparation; includes upper and response construction and repeated exact replays.')
    if temple_offset is not None:
        r['discovery_scope'] = 'Old P and reference factors are reused; upper and all response directions are generated on the connected Hamiltonian. The response targets a fixed excitation count, followed by exact Temple replay. No inherited response prefix, joined-Q factorization, or general convergence guarantee.'
    if proof_key=='complement_atoms':
        r['reference_gap_family']=base['complement_atoms']['kind']
        r['discovery_scope']='Old P and complete reference atom/DD proof are reused. The physical reference witness seeds a new upper refinement; every response direction is discovered on the perturbed Hamiltonian. Exact norm-shift, Schur and final energy gates remain mandatory. No dense Q factor or scalable frontier is claimed.'
    (out / 'certificate.json').write_text(json.dumps(final, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(r, indent=2) + '\n')
    print(json.dumps({k:r[k] for k in ('width_float','response_dimension','current_q_component_dimensions','largest_factor_dimension','elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--strength', default='1/1000')
    parser.add_argument('--upper-steps', type=int, default=24)
    parser.add_argument('--norm-squares', action='store_true')
    parser.add_argument('--temple-offset')
    args = parser.parse_args()
    transfer(args.source, args.output, F(args.strength), args.upper_steps, norm_squares=args.norm_squares,
             temple_offset=F(args.temple_offset) if args.temple_offset is not None else None)
