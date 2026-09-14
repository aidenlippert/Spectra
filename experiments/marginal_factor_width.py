"""Exact dual obstruction to configuration-supported factor-width positivity.

The support size is in the fixed determinant coordinate basis. It is not
fermion body degree, orbital locality, or a lower bound for every proof family.
"""
from fractions import Fraction as F
import json
from pathlib import Path

from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_sparse_response import parse_basis
from experiments.marginal_spin_reduction import SpinZeroOracle


def pairing(oracle, retained, vector, width, threshold):
    """Trace((QHQ-threshold)B)/Tr(B), for the implicit k-positive dual B."""
    norm = sum(x*x for x in vector.values())
    value = F(0)
    for s, x in vector.items():
        for t, coefficient in oracle.action(s).items():
            if t not in retained and t in vector:
                value += x*vector[t]*(coefficient if s == t else -abs(coefficient)/F(width-1))
    return value/norm-threshold


def replay(certificate):
    if certificate.get('kind') != 'spin_factor_width_obstruction_v1':
        raise ValueError('Unsupported factor-width obstruction')
    width = certificate.get('factor_width')
    if type(width) is not int or not 2 <= width <= 64:
        raise ValueError('Factor width must be an integer from two to64')
    oracle = SpinZeroOracle(certificate)
    retained = oracle.retained(certificate.get('retained_states'))
    vector = parse_basis(oracle, retained, [certificate.get('dual_weights')])[0]
    if any(x <= 0 for x in vector.values()):
        raise ValueError('Strictly positive dual weights required')
    threshold = F(certificate['target_lower'])
    value = pairing(oracle, retained, vector, width, threshold)
    if value >= 0:
        raise ValueError('Dual pairing does not obstruct the requested factor width')
    return {'factor_width_excluded': width, 'minimum_possible_maximum_support': width+1,
        'target_lower': rational_text(threshold), 'normalized_dual_pairing': rational_text(value),
        'normalized_dual_pairing_float': float(value), 'witness_support': len(vector),
        'unique_action_states': len(oracle.cache), 'referenced_determinants': oracle.referenced_state_count(),
        'scope': 'Exact negative trace pairing with a dual matrix whose every principal submatrix of order at most k is PSD. Excludes every sum of PSD pieces supported on at most k determinant coordinates for QHQ-target. No complete Q enumeration in replay. This is not an obstruction to physical positivity, arbitrary coordinate changes, low fermion body degree, or other certificate representations.'}


def export(source, output, width=3):
    """Small explicit block eigenvectors propose weights; replay is independent."""
    import numpy as np
    from experiments.marginal_h6_complement import checked_blocks
    source, out = Path(source), Path(output)
    if type(width) is not int or not 2 <= width <= 64:
        raise ValueError('Factor width must be an integer from two to64')
    if out.exists():
        raise ValueError('Preserve previous factor-width export')
    c = json.loads(source.read_text())['spin_symmetric_certificate']
    oracle = SpinZeroOracle(c)
    retained = oracle.retained(c['retained_states'])
    matrices = checked_blocks(oracle, c['retained_states'], c['blocks'])
    threshold = F(c['complement_lower'])
    history = []
    accepted = None
    for block, a in zip(c['blocks'], matrices):
        n = len(a)
        matrix = -np.abs(np.array(a, dtype=float))/(width-1)
        np.fill_diagonal(matrix, np.diag(np.array(a, dtype=float))-float(threshold))
        values, vectors = np.linalg.eigh(matrix)
        history.append({'block_dimension': n, 'numerical_smallest_eigenvalue': float(values[0])})
        if values[0] >= 0:
            continue
        weights = np.abs(vectors[:, 0])
        order = sorted(range(n), key=lambda i: (-weights[i], block['states'][i]))
        for count in sorted(set(min(n,k) for k in (4,8,16,24,32,48,64,96,128,160,200,256))):
            chosen = sorted(order[:count], key=lambda i: block['states'][i])
            recipe = {'states': [block['states'][i] for i in chosen],
                      'amplitudes': [int(round(weights[i]*10**10)) for i in chosen]}
            trial = {key: c[key] for key in ('modes','particles','hamiltonian','retained_states')}
            trial.update(kind='spin_factor_width_obstruction_v1', factor_width=width,
                         target_lower=rational_text(threshold), dual_weights=recipe)
            try:
                receipt = replay(trial)
            except ValueError:
                continue
            accepted = trial, receipt
            break
        if accepted is not None:
            break
    if accepted is None:
        raise ValueError('This restricted dual family found no exact obstruction')
    certificate, receipt = accepted
    receipt.update(source=str(source), discovery_unique_action_states=len(oracle.cache),
        discovery_scope='Uses explicit reference blocks and their comparison eigenvectors to propose a dual witness; final replay uses only the exported witness source states. Failure of this proposer would not prove membership in the factor-width cone.')
    out.mkdir(parents=True)
    for name, value in (('certificate.json',certificate),('receipt.json',receipt),('proposal.json',history)):
        (out/name).write_text(json.dumps(value,indent=2)+'\n')
    print(json.dumps({key:receipt[key] for key in ('factor_width_excluded','normalized_dual_pairing_float','witness_support','unique_action_states')}),flush=True)


def replay_transfer(certificate):
    """Exclude an FW-k reference gap with a saturated scalar norm shift."""
    from experiments.marginal_operator_norm import replay as norm_replay, replay_q
    from experiments.marginal_spin_temple import replay as interval_replay
    from experiments.marginal_symbolic import add, scale
    flat=certificate.get('kind')=='spin_flat_factor_transfer_obstruction_v1'
    if not flat and certificate.get('kind') != 'spin_factor_width_transfer_obstruction_v1':
        raise ValueError('Unsupported factor-width transfer obstruction')
    reference_certificate = certificate.get('reference_obstruction')
    if flat:
        from experiments.marginal_clique_gap import replay_flat
        reference_receipt=replay_flat(reference_certificate)
        excluded=reference_receipt['maximum_support_excluded']
    else:
        reference_receipt=replay(reference_certificate)
        excluded=reference_receipt['factor_width_excluded']
    reference = SpinZeroOracle(reference_certificate)
    interval = certificate.get('current_interval')
    current_data = interval.get('spin_symmetric_certificate') if type(interval) is dict else None
    if (type(current_data) is not dict or current_data.get('modes') != reference.modes
            or current_data.get('particles') != reference.particles
            or current_data.get('retained_states') != reference_certificate['retained_states']):
        raise ValueError('Current excitation proof and reference obstruction must use the same sector and P')
    current = SpinZeroOracle(current_data)
    norm_certificate = certificate.get('saturated_norm')
    if (type(norm_certificate) is not dict or norm_certificate.get('modes') != reference.modes
            or norm_certificate.get('particles') != reference.particles):
        raise ValueError('Saturated norm must use the same sector')
    delta = SpinZeroOracle({'modes':reference.modes, 'particles':reference.particles,
                           'hamiltonian':norm_certificate['perturbation']})
    if delta.h != add(current.h, scale(reference.h,-1)):
        raise ValueError('Norm perturbation must equal the actual current-reference difference')
    compressed = norm_certificate.get('kind') == 'saturated_q_operator_norm_v1'
    if compressed and norm_certificate.get('retained_states') != reference_certificate['retained_states']:
        raise ValueError('Compressed norm must use the same retained P')
    norm_receipt = replay_q(norm_certificate) if compressed else norm_replay(norm_certificate)
    norm = F(norm_receipt['operator_norm'])
    # Saturation in Sz=0 prevents an improvement merely by restricting the
    # full-sector norm. The optional Q receipt also checks saturation inside Q.
    if abs(delta.upper(norm_certificate['saturating_witness'])) != norm:
        raise ValueError('Norm witness must saturate in the same spin sector')
    interval_receipt = interval_replay(interval)
    current_lower = F(interval_receipt['spin_symmetric_lower'])
    ceiling = F(reference_receipt['target_lower'])+F(reference_receipt['normalized_dual_pairing'])
    shifted = ceiling-norm
    if shifted > current_lower:
        raise ValueError('The factor-width ceiling does not obstruct every valid upper witness')
    result = {'factor_width_excluded': excluded,
        'reference_threshold_ceiling': rational_text(ceiling),
        'exact_full_spin_sector_perturbation_norm': rational_text(norm),
        'best_possible_shifted_threshold_ceiling': rational_text(shifted),
        'current_spin_symmetric_ground_lower': rational_text(current_lower),
        'threshold_deficit_below_ground_lower': rational_text(current_lower-shifted),
        'threshold_deficit_float': float(current_lower-shifted),
        'reference_dual_witness_support': reference_receipt['witness_support'],
        'current_interval_unique_source_states': interval_receipt['unique_determinant_sources'],
        'scope': 'For this reference H0, fixed P, determinant basis, and current Hs, every factor-width-k Q-gap certificate shifted by any valid full-Sz=0 perturbation norm falls at or below a certified current ground lower. It therefore cannot exceed any valid upper witness. Does not exclude Q-compressed norm improvements, a changed P/reference/basis, wider factors, or other positivity proofs. Replays a complete current interval in addition to the compact dual witness.'}
    if compressed:
        result['exact_q_compressed_perturbation_norm'] = rational_text(norm)
        result['scope'] = 'For this reference H0, fixed P, determinant basis, and current Hs, every factor-width-k Q-gap certificate shifted by any valid Q-compressed perturbation norm falls at or below a certified current ground lower. Thus this scalar norm-shift pipeline cannot exceed any valid upper witness, even with the optimal compressed norm. Does not exclude a changed P/reference/basis, wider factors, directional perturbation bounds, or other positivity proofs. Replays the current interval, dual witness, and exact Q-norm saturation.'
    if flat:
        result['equal_magnitude_support_excluded']=result.pop('factor_width_excluded')
        result['reference_obstruction_source_states']=reference_receipt['unique_action_states']
        result['scope']='For this H0, P, basis, and Hs, every Q-gap proof made from equal-magnitude rank-one squares of the stated maximum support, followed by the certified saturated scalar norm shift, falls below a certified current ground lower. With the Q-norm receipt this includes the optimal compressed norm. Arbitrary-amplitude factors, larger supports, coordinate changes, and directional bounds remain outside the excluded family.'
    return result


def replay_separator(certificate):
    """Check a (k+1)-coordinate positive atom violated by the k-positive dual."""
    if certificate.get('kind') != 'spin_factor_width_separator_v1':
        raise ValueError('Unsupported factor-width separator')
    parent = certificate.get('dual_obstruction')
    parent_receipt = replay(parent)
    oracle = SpinZeroOracle(parent)
    retained = oracle.retained(parent['retained_states'])
    weights = parse_basis(oracle, retained, [parent['dual_weights']])[0]
    vector = parse_basis(oracle, retained, [certificate.get('separator')])[0]
    width = parent['factor_width']
    if len(vector) != width+1 or not set(vector).issubset(weights):
        raise ValueError('Separator must use k+1 coordinates in the dual support')
    value = sum((x*weights[s])**2 for s,x in vector.items())
    for s,x in vector.items():
        for t,coefficient in oracle.action(s).items():
            if t != s and t in vector:
                sign = (coefficient>0)-(coefficient<0)
                value -= F(sign,width-1)*x*weights[s]*vector[t]*weights[t]
    norm = sum(x*x for x in vector.values())
    if value >= 0:
        raise ValueError('Positive atom does not separate the exported dual witness')
    normalization = norm*sum(x*x for x in weights.values())
    return {'atom_support':len(vector), 'excluded_dual_factor_width':width,
        'normalized_dual_atom_pairing':rational_text(value/normalization),
        'normalized_dual_atom_pairing_float':float(value/normalization),
        'separator_check_source_states':len(oracle.cache),
        'parent_obstruction_source_states':parent_receipt['unique_action_states'],
        'scope':'The rank-one PSD atom zz^T on k+1 determinant coordinates has negative exact pairing with the exported k-positive dual; both matrices are normalized to unit trace in the reported pairing. Adding this constraint excludes that dual witness; it does not prove factor-width-(k+1) membership or sufficiency. Complete replay also verifies the parent obstruction.'}


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify')
    parser.add_argument('--verify-transfer')
    parser.add_argument('--verify-separator')
    parser.add_argument('--source')
    parser.add_argument('--output')
    parser.add_argument('--width',type=int,default=3)
    args=parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())),indent=2))
    elif args.verify_transfer:
        print(json.dumps(replay_transfer(json.loads(Path(args.verify_transfer).read_text())),indent=2))
    elif args.verify_separator:
        print(json.dumps(replay_separator(json.loads(Path(args.verify_separator).read_text())),indent=2))
    elif args.source and args.output:
        export(args.source,args.output,args.width)
    else:
        parser.error('Specify --verify or --source and --output')
