"""Exact SU(2) ground-sector reduction with a checked rounding perturbation."""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_h6_complement import verify_complement
from experiments.marginal_implicit_certificate import rational_text
from experiments.marginal_schur_transfer import ldl_pivots
from experiments.marginal_sparse_response import parse_basis, prepare, response_matrix
from experiments.marginal_symbolic import add, adj, mono, product, scale, canonical, encode


def commutator(a, b):
    return add(product(a, b), scale(product(b, a), -1))


def spin_operators(modes):
    if type(modes) is not int or modes < 2 or modes % 2:
        raise ValueError('Paired alpha/beta spin orbitals required')
    plus = add(*(mono(((1, 2*i), (0, 2*i+1))) for i in range(modes // 2)))
    z = add(*(mono(((1, i), (0, i)), F(1 if i % 2 == 0 else -1, 2)) for i in range(modes)))
    return plus, canonical(adj(plus)), z


def symmetrize(h, modes):
    plus, minus, z = spin_operators(modes)
    def casimir(a):
        return add(commutator(z, commutator(z, a)),
                   scale(add(commutator(plus, commutator(minus, a)),
                             commutator(minus, commutator(plus, a))), F(1, 2)))
    # Operator spin ranks 0,1,2 have Casimir eigenvalues 0,2,6.
    # The proposer uses this projection; replay independently checks commutators.
    first = casimir(h)
    return canonical(add(h, scale(first, F(-2, 3)), scale(casimir(first), F(1, 12))))


class SpinZeroOracle(DeterminantOracle):
    def __init__(self, certificate):
        super().__init__(certificate)
        if self.particles % 2 or self.modes % 2:
            raise ValueError('Spin-zero reduction requires even modes and even particle number')
        self.alpha_mask = sum(1 << i for i in range(0, self.modes, 2))
        plus, _, z = spin_operators(self.modes)
        if commutator(self.h, plus) or commutator(self.h, z):
            raise ValueError('Exact spin commutation is required; approximate symmetry is insufficient')

    def valid_state(self, state):
        return super().valid_state(state) and (state & self.alpha_mask).bit_count() == self.particles // 2

    @property
    def sector_dimension(self):
        return math.comb(self.modes // 2, self.particles // 2)**2

    def cover(self, *args, **kwargs):
        raise ValueError('Use explicit spin-sector block coverage, not the fixed-N occupation tree')


def gap_coverage(oracle, p, gamma, recipe):
    if 'complement_atoms' not in recipe:
        return verify_complement(oracle,p,gamma,recipe.get('blocks'))
    if 'blocks' in recipe:
        raise ValueError('Choose one direct complement proof family')
    atoms=recipe['complement_atoms']
    if type(atoms) is not dict or atoms.get('kind') not in ('spin_signed_atom_complement_v1','spin_rational_atom_complement_v1','spin_charge_product_complement_v1','valence_coherent_tree_v1','valence_charge_spin_v1','joint_polynomial_metric_gap_v1','joint_polynomial_monotone_transfer_v1'):
        raise ValueError('Explicit supported complement atom recipe required')
    original_kind=atoms['kind']
    if original_kind in ('valence_coherent_tree_v1','valence_charge_spin_v1','joint_polynomial_metric_gap_v1','joint_polynomial_monotone_transfer_v1'):
        if original_kind=='valence_coherent_tree_v1':
            from experiments.marginal_coherent_tree import replay as replay_gap
            allowed={'kind','metric_rule','tree','max_nodes','diagonal_charge_lambda','source_polynomial_bound','joint_line_bound','joint_simplex_bound'}
        elif original_kind=='valence_charge_spin_v1':
            from experiments.marginal_charge_spin import replay as replay_gap
            allowed={'kind','metric_rule','max_charge_patterns'}
        elif original_kind=='joint_polynomial_monotone_transfer_v1':
            from experiments.marginal_monotone_transfer import replay as replay_gap
            allowed={'kind','reference_certificate'}
        else:
            from experiments.marginal_polynomial_metric import replay as replay_gap
            allowed={'kind','polynomial_metric','weight_proof','numerator_proof'}
        if set(atoms)-allowed:
            raise ValueError('Charge proof recipe contains unsupported or hidden bound data')
        retained=oracle.retained(p)
        if (oracle.particles!=oracle.modes//2 or oracle.modes%4
                or len(retained)!=math.comb(oracle.modes//2,oracle.particles//2)
                or any(any(((s>>(2*i))&3) not in (1,2) for i in range(oracle.modes//2)) for s in retained)):
            raise ValueError('Charge proof requires the complete valence reference')
        bound=dict(atoms,modes=oracle.modes,particles=oracle.particles,
                   hamiltonian=encode(oracle.h),target_lower=str(F(gamma)))
        result=replay_gap(bound)
        result['proof_family']=original_kind
        return result
    if original_kind=='spin_charge_product_complement_v1':
        from experiments.marginal_charge_product import expand_recipe
        atoms=expand_recipe(atoms,oracle.modes)
    from experiments.marginal_h6_complement import checked_blocks
    from experiments.marginal_signed_atoms import verify_block
    blocks=atoms.get('blocks')
    matrices=checked_blocks(oracle,p,blocks,max_block_dimension=384)
    rational=atoms['kind']=='spin_rational_atom_complement_v1'
    rows=[verify_block(a,b,gamma,rational) for a,b in zip(matrices,blocks)]
    return {'complement_lower':rational_text(F(gamma)),'blocks':rows,
            'proof_family':original_kind,
            'scope':'Positive local atoms and exact complete weighted DD residual, rebuilt from this Hamiltonian and retained space; no dense complement factor.'}


def spin_gap(oracle, p, reduced):
    gamma = F(reduced['complement_lower'])
    if 'complement_reference' not in reduced:
        return gap_coverage(oracle, p, gamma, reduced), None
    if 'complement_atoms' in reduced:
        raise ValueError('Choose direct atoms or a transferred complement reference')
    reference = reduced['complement_reference']
    if type(reference) is not dict:
        raise ValueError('Explicit complement reference required')
    other = SpinZeroOracle({'modes': oracle.modes, 'particles': oracle.particles,
                            'hamiltonian': reference['hamiltonian']})
    threshold = F(reference['complement_lower'])
    coverage = gap_coverage(other, p, threshold, reference)
    difference = add(oracle.h, scale(other.h, -1))
    norm_receipt = None
    if 'norm_certificate' in reference:
        from experiments.marginal_operator_norm import verify_norm
        error, norm_receipt = verify_norm(difference, oracle.modes, reference['norm_certificate'])
    else:
        error = sum(abs(x) for x in difference.values())
    if gamma > threshold - error:
        raise ValueError('Requested complement threshold exceeds the reference norm-shift bound')
    result = {'complement_lower': rational_text(gamma), 'reference_coverage': coverage,
            'reference_perturbation_norm_bound': rational_text(error),
            'shifted_reference_threshold': rational_text(threshold - error)}
    if norm_receipt is not None:
        result['norm_square_receipt'] = norm_receipt
    return result, other


def replay(certificate, *, audit=None):
    if certificate.get('kind') != 'spin_reduced_perturbation_interval_v1':
        raise ValueError('Unsupported spin reduction certificate')
    original = DeterminantOracle(certificate)
    reduced = certificate.get('spin_symmetric_certificate')
    if (type(reduced) is not dict or reduced.get('modes') != original.modes
            or reduced.get('particles') != original.particles):
        raise ValueError('Spin reference must describe the same physical sector')
    oracle = SpinZeroOracle(reduced)
    p = reduced.get('retained_states')
    retained = oracle.retained(p)
    gamma, lower = F(reduced['complement_lower']), F(reduced['lower'])
    if gamma <= lower:
        raise ValueError('Positive spin-sector complement denominator required')
    coverage, reference_oracle = spin_gap(oracle, p, reduced)
    basis = parse_basis(oracle, retained, reduced.get('response_basis'), 32)
    if ldl_pivots(response_matrix(prepare(oracle, p, gamma, basis), lower)) is None:
        raise ValueError('Spin-sector Schur positivity failed')
    error = sum(abs(x) for x in add(original.h, scale(oracle.h, -1)).values())
    # Each CAR monomial has norm <=1. Exact coefficient l1 gives a full-space error bound.
    lower -= error
    witness = reduced.get('independent_upper')
    oracle.upper(witness)  # Enforce the advertised spin-zero support.
    upper = original.upper(witness)
    if lower > upper:
        raise ValueError('Inconsistent original-Hamiltonian interval')
    oracles = [original, oracle] + ([reference_oracle] if reference_oracle is not None else [])
    source_states = set().union(*(set(o.cache) for o in oracles))
    referenced = source_states | {s for o in oracles for row in o.cache.values() for s in row}
    if audit is not None:
        audit['sources'].update(source_states)
        audit['referenced'].update(referenced)
    result = {'lower': rational_text(lower), 'upper': rational_text(upper),
            'width': rational_text(upper - lower), 'width_float': float(upper - lower),
            'spin_symmetry_error_bound': rational_text(error),
            'spin_symmetric_lower': reduced['lower'], 'retained_dimension': len(p),
            'response_dimension': len(basis), 'spin_sector_dimension': oracle.sector_dimension,
            'full_sector_dimension': original.sector_dimension, 'complement_coverage': coverage,
            'unique_determinant_sources': len(source_states), 'referenced_determinants': len(referenced),
            'spin_symmetric_action_states': len(oracle.cache), 'original_action_states': len(original.cache),
            'scope': 'Global original-Hamiltonian interval from exact SU(2) spin-zero reduction of a nearby Hamiltonian, complete spin-sector Q factor coverage, and a coefficient-norm perturbation bound. Ground-energy reduction does not assert symmetry of P or of the original projected Q. All spin-zero configurations are explicit; no polynomial scaling claim.'}
    if reference_oracle is not None:
        result['reference_action_states'] = len(reference_oracle.cache)
        result['scope'] = 'Global original-Hamiltonian interval with exact SU(2) reduction, a complete reference Q factor proof shifted by a recomputed operator-norm bound, and residual response of the actual spin-symmetric Hamiltonian. Original-H transfer is checked separately. Current Q may connect the reference blocks; no factorization of the joined block or asymptotic scaling claim.'
    return result


def reduce_certificate(path, output):
    source, out = Path(path), Path(output)
    if out.exists():
        raise ValueError('Preserve previous spin reduction export')
    c = json.loads(source.read_text())
    if c.get('kind') != 'factor_response_schur_v1':
        raise ValueError('Expected a factor-response source')
    started = time.monotonic()
    original = DeterminantOracle(c)
    h = symmetrize(original.h, original.modes)
    reduced = {key: c[key] for key in ('modes', 'particles', 'retained_states', 'complement_lower', 'lower', 'response_basis', 'independent_upper')}
    reduced['hamiltonian'] = encode(h)
    oracle = SpinZeroOracle(reduced)
    # Reuse factors only for entire blocks lying in this exact conserved sector.
    blocks = []
    multiplier = math.lcm(*(x.denominator for x in h.values()), F(c['complement_lower']).denominator)
    for item in c['blocks']:
        included = [oracle.valid_state(s) for s in item['states']]
        if not any(included):
            continue
        if not all(included):
            raise ValueError('A source factor block straddles the chosen spin sector')
        f = item['factor']
        ratio = multiplier // math.gcd(multiplier, f['scale'])
        blocks.append({'states': item['states'], 'factor': {
            'scale': f['scale'] * ratio, 'diagonal': [x * ratio for x in f['diagonal']],
            'lower': [[x * ratio for x in row] for row in f['lower']]}})
    reduced['blocks'] = blocks
    result = {key: c[key] for key in ('modes', 'particles', 'hamiltonian')}
    result.update(kind='spin_reduced_perturbation_interval_v1', spin_symmetric_certificate=reduced)
    error = sum(abs(x) for x in add(original.h, scale(h, -1)).values())
    for safety in (0, 1, 2, 10, 100):
        reduced['lower'] = rational_text(F(c['lower']) - safety * max(error, F(1, 10**12)))
        try:
            receipt = replay(result)
        except ValueError:
            continue
        break
    else:
        raise ValueError('Inherited spin-sector factors/response failed exact replay')
    receipt.update(source=str(source), elapsed_seconds=time.monotonic() - started,
                   discovery_scope='Inherits an existing full-Q factor-response certificate; only final replay has reduced state coverage. End-to-end reduced construction is not claimed.')
    out.mkdir(parents=True)
    (out / 'certificate.json').write_text(json.dumps(result, indent=2) + '\n')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({key: receipt[key] for key in ('width_float', 'spin_symmetry_error_bound', 'spin_sector_dimension', 'unique_determinant_sources', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify')
    parser.add_argument('--reduce')
    parser.add_argument('--output')
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(replay(json.loads(Path(args.verify).read_text())), indent=2))
    elif args.reduce and args.output:
        reduce_certificate(args.reduce, args.output)
    else:
        parser.error('Specify --verify or --reduce and --output')
