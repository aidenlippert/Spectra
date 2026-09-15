"""Exact Casimir, common-kernel energy, and charged spin-separator diagnostics."""
from fractions import Fraction as F
from math import comb, lcm
from pathlib import Path
import hashlib
import json
import sys
import time

from experiments.marginal_symbolic import add, canonical, mono, product, scale
from experiments.marginal_hunt_car import adj
from research.molecular_collective_20260913.core import extract, density, tail_replay, retained_polynomial, digest
from research.certificate_scaling.streaming_reference_upper import upper
from research.certificate_scaling.commutator_dual_witness import moment_decode, evaluate, seed
from research.response_consistency_20260913.closure import lie_closure
from research.response_consistency_20260913 import separator
from research.positive_response_20260913.coercivity import spatial_inputs, check as gap_check

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/positive_response_20260913'


def commutator(a, b):
    return add(product(a, b), scale(product(b, a), -1))


def spin_operators(s):
    z = add(*(mono(((1, 2*i+sigma), (0, 2*i+sigma)), F(1-2*sigma, 2))
              for i in range(s) for sigma in range(2)))
    plus = add(*(mono(((1, 2*i), (0, 2*i+1))) for i in range(s)))
    minus = canonical(adj(plus))
    squared = add(product(z, z), scale(add(product(plus, minus), product(minus, plus)), F(1, 2)))
    return z, plus, minus, squared


def casimir_identity(s):
    start = time.monotonic(); number = add(*(density(i, i) for i in range(s)))
    _, _, _, spin_squared = spin_operators(s)
    terms = [add(density(i, j), scale(number, -F(1, s))) if i == j else density(i, j)
             for i in range(s) for j in range(s)]
    orbital = add(*(product(canonical(adj(op)), op) for op in terms))
    rhs = add(scale(number, s+2), scale(product(number, number), -F(1, 2)-F(1, s)), scale(spin_squared, -2))
    if orbital != rhs:
        raise AssertionError('Orbital/spin Casimir identity failed in exact CAR arithmetic')
    return {'spatial_orbitals': s, 'all_particle_sectors_CAR_identity': True,
        'degree_four_polynomial_terms': len(orbital), 'many_body_states_enumerated': 0,
        'wall_seconds': time.monotonic()-start}


def spatial_one_body(p):
    s = p['spatial']; one = p['one']
    matrix = [[one.get(((1, 2*i), (0, 2*j)), F(0)) for j in range(s)] for i in range(s)]
    constant = one.get((), F(0))
    rebuilt = add(mono((), constant), *(scale(density(i, j), matrix[i][j]) for i in range(s) for j in range(s)))
    if rebuilt != one or matrix != list(map(list, zip(*matrix))):
        raise ValueError('Exact spin-independent real symmetric one-body term required')
    return constant, matrix


def complement_bounds(constant, matrix, bulk_gap, tail_lower, target):
    s = len(matrix); remaining_particles = s-2
    if remaining_particles % 2:
        raise ValueError('Even half filling required for this paired diagonal lower bound')
    rows = []
    for orbital in range(s):
        others = [i for i in range(s) if i != orbital]
        rho = max(sum(abs(matrix[i][j]) for j in others if j != i) for i in others)
        diagonal = sorted(matrix[i][i] for i in others)
        one_lower = constant+2*matrix[orbital][orbital]+2*sum(diagonal[:remaining_particles//2])-remaining_particles*rho
        lower = one_lower+bulk_gap+tail_lower; delta = lower-target
        without_bulk = one_lower+tail_lower-target
        rows.append({'doubly_occupied_spatial_orbital': orbital,
            'one_body_lower_Ha': str(one_lower), 'remaining_spatial_offdiagonal_row_bound_Ha': str(rho),
            'full_complement_lower_Ha': str(lower), 'shifted_D_lower_Ha': str(delta),
            'shifted_D_lower_float_Ha': float(delta), 'strict_gap_certified': delta > 0,
            'shifted_D_lower_without_commutator_gap_Ha': str(without_bulk),
            'shifted_D_lower_without_commutator_gap_float_Ha': float(without_bulk),
            'one_body_plus_R_nonnegative_already_suffices': without_bulk > 0,
            'eliminated_sector_dimension_formula': f'binomial({2*s-2},{s-2})',
            'eliminated_sector_dimension': comb(2*s-2, s-2), 'many_body_states_enumerated': 0})
    return rows


def sector_gap_check(data, tail, cert):
    """A standalone molecular D gap using R>=0; no Lie reconstruction needed."""
    start = time.monotonic()
    keys = {'kind', 'fixture_sha256', 'tail_sha256', 'doubly_occupied_spatial_orbital', 'target_Ha', 'delta_Ha'}
    if set(cert) != keys or cert['kind'] != 'molecular_double_occupancy_gap_v1':
        raise ValueError('Unknown molecular sector-gap certificate')
    if cert['fixture_sha256'] != digest(data) or cert['tail_sha256'] != digest(tail):
        raise ValueError('Sector-gap input binding failed')
    if data['modes'] != 12 or data['particles'] != 6:
        raise ValueError('This molecular sector-gap rule is capped at H6')
    orbital = cert['doubly_occupied_spatial_orbital']
    if type(orbital) is not int or not 0 <= orbital < 6:
        raise ValueError('Invalid eliminated double-occupancy subspace')
    if any(type(cert[k]) is not str for k in ('target_Ha', 'delta_Ha')):
        raise ValueError('Exact rational target and gap required')
    target = F(cert['target_Ha']); delta = F(cert['delta_Ha'])
    if delta <= 0:
        raise ValueError('Strictly positive eliminated-sector gap required')
    p = extract(data, tail['center_number']); constant, matrix = spatial_one_body(p)
    tail_receipt = tail_replay(data, tail); ell = F(tail_receipt['lower_operator_shift_Ha'])
    row = complement_bounds(constant, matrix, F(0), ell, target)[orbital]
    if F(row['shifted_D_lower_Ha']) < delta:
        raise ValueError('Claimed eliminated-sector gap exceeds the exact lower bound')
    return row | {'accepted_delta_Ha': str(delta), 'accepted_delta_float_Ha': float(delta),
        'target_Ha': str(target), 'tail_replay': tail_receipt,
        'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'many_body_inverse_entries': 0, 'lie_bracket_nodes_required': 0,
        'replay_seconds': time.monotonic()-start,
        'scope': f'Q(H-b)Q>=delta Q for Q=n_({orbital},up)n_({orbital},down). This does not certify H>=b or a response.'}


def molecular_case(name, rank, gap_cert=None):
    start = time.monotonic(); src = ROOT/'results/molecular_collective_20260913/campaign'/name
    data = json.loads((src/'fixture.json').read_text()); tail = json.loads((src/f'rank_{rank}/tail.json').read_text())
    p = extract(data, tail['center_number']); s = p['spatial']
    matrices, _, _ = spatial_inputs(data, tail)
    before = time.monotonic(); closure = lie_closure(matrices); closure['wall_seconds'] = time.monotonic()-before
    if not closure['full_sl_closure_proved']:
        raise ValueError('Common kernel theorem requires full traceless closure')
    casimir = casimir_identity(s); constant, matrix = spatial_one_body(p)
    tail_receipt = tail_replay(data, tail); tail_lower = F(tail_receipt['lower_operator_shift_Ha'])
    kernel_energy = constant+sum(matrix[i][i] for i in range(s))
    all_up = sum(1 << (2*i) for i in range(s))
    if separator.act(retained_polynomial(p, tail), {all_up: F(1)}) != {all_up: kernel_energy}:
        raise AssertionError('Independent determinant action disagrees with the kernel energy')
    upper_path = ROOT/'results/certificate_scaling/active_space_ladder_references_aligned'/name/'upper.json'
    upper_bytes = upper_path.read_bytes(); witness = json.loads(upper_bytes)
    U, upper_receipt = upper(data, witness['independent_upper'])
    if U != F(witness['upper']):
        raise AssertionError('Frozen physical upper changed')
    exclusion = kernel_energy+tail_lower-U
    if exclusion <= 0:
        raise ValueError('This diagnostic did not exclude the common kernel from the ground space')
    result = {'molecule': name, 'fixture_sha256': digest(data), 'tail_sha256': digest(tail),
        'upper_witness_file_sha256': hashlib.sha256(upper_bytes).hexdigest(),
        'retained_patterns': len(matrices), 'lie_closure': closure, 'casimir_identity': casimir,
        'maximum_spin': str(F(s, 2)), 'common_kernel_dimension': s+1,
        'retained_kernel_energy_Ha': str(kernel_energy), 'retained_kernel_energy_float_Ha': float(kernel_energy),
        'original_H_kernel_lower_Ha': str(kernel_energy+tail_lower),
        'physical_upper_Ha': str(U), 'physical_upper_float_Ha': float(U),
        'kernel_energy_above_physical_upper_Ha': str(exclusion),
        'kernel_energy_above_physical_upper_float_Ha': float(exclusion),
        'common_kernel_excluded_from_ground_space': True, 'tail_replay': tail_receipt,
        'upper_replay': upper_receipt, 'new_determinants_generated': 1,
        'inherited_upper_amplitudes_read': upper_receipt['witness_states'],
        'upper_discovery_cost': 'Inherited FCI-based witness, 200 H6 / 1000 H8 amplitudes; no discovery performed in this pass.'}
    if gap_cert is not None:
        gap = gap_check(data, tail, gap_cert); bulk_gap = F(gap['double_occupancy_bulk_gap_Ha'])
        target = U-F(1, 625)  # Exactly 1.6 mHa below the frozen physical upper.
        result['gap_replay'] = gap
        result['complement_gap_target_Ha'] = str(target)
        result['complement_gap_target_float_Ha'] = float(target)
        result['complement_gap_bounds'] = complement_bounds(constant, matrix, bulk_gap, tail_lower, target)
    result['wall_seconds'] = time.monotonic()-start
    result['scope'] = 'Excludes a common-zero model of these retained squares at low molecular energy. Does not exclude shifted or frustrated positive decompositions or other response constructions.'
    return result


def spin_half(C, s):
    z, plus, minus, _ = spin_operators(s)
    def casimir(op):
        return add(commutator(z, commutator(z, op)),
            scale(add(commutator(plus, commutator(minus, op)), commutator(minus, commutator(plus, op))), F(1, 2)))
    projected = scale(add(scale(C, F(15, 4)), scale(casimir(C), -1)), F(1, 3))
    if add(commutator(z, C), scale(C, F(1, 2))):
        raise ValueError('Frozen source must be a spin component with m=-1/2')
    if not projected or add(casimir(projected), scale(projected, -F(3, 4))):
        raise AssertionError('Spin-one-half Casimir eigenvalue failed')
    if commutator(minus, projected) or add(commutator(z, projected), scale(projected, F(1, 2))):
        raise AssertionError('Spin-doublet lowest-component covariance failed')
    return projected


def projected_separator(source_bytes, parent_bytes):
    source = json.loads(source_bytes); C, _ = separator.operator(source); projected = spin_half(C, 6)
    den = lcm(*(c.denominator for c in projected.values())); words = sorted(projected)
    cert = {'kind': 'three_removal_separator_v1', 'modes': 12,
        'parent_witness_sha256': hashlib.sha256(parent_bytes).hexdigest(),
        'source_separator_sha256': hashlib.sha256(source_bytes).hexdigest(),
        'projection': '(15*C/4-J_ad_squared(C))/3', 'support_spatial': sorted({i//2 for w in words for _, i in w}),
        'triples': [[i for _, i in w] for w in words], 'coefficients': [int(projected[w]*den) for w in words],
        'denominator': den}
    _, P = separator.operator(cert); y = moment_decode(json.loads(parent_bytes)['moments'], 12)
    cert['exact_moment'] = str(evaluate(P, y))
    return cert


def check_projected(cert, source_bytes, parent_bytes):
    start = time.monotonic(); expected = projected_separator(source_bytes, parent_bytes)
    if cert != expected:
        raise ValueError('Certificate is not the exact hash-bound spin projection')
    receipt = separator.check(cert, parent_bytes); C, P = separator.operator(cert)
    if receipt['sparse_zero_terms_used']:
        raise ValueError('All projected-separator moments must be explicitly present')
    z, plus, minus, _ = spin_operators(6); raised = commutator(plus, C)
    # Restore the omitted spin-one-half partner and test the spin-scalar sum too.
    from research.joint_patterns_20260913.core import anticommutator
    P_partner = anticommutator(raised, raised); P_scalar = add(P, P_partner)
    if commutator(z, P_scalar) or commutator(plus, P_scalar) or commutator(minus, P_scalar):
        raise AssertionError('Completed doublet did not form a spin scalar')
    y = moment_decode(json.loads(parent_bytes)['moments'], 12)
    if any(w not in separator.domain(12) or w not in y for w in P_scalar):
        raise ValueError('Spin completion requires unspecified moments')
    scalar_moment = evaluate(P_scalar, y)
    if scalar_moment >= 0:
        raise ValueError('Spin-scalar completion did not reject the frozen functional')
    vector = {sum(1 << (2*i) for i in range(6)): F(1)}; controls = []
    for step in range(7):
        if not vector or separator.act(C, vector) or separator.act(canonical(adj(C)), vector):
            raise AssertionError('Projected operator failed to annihilate the maximum-spin multiplet')
        if separator.act(raised, vector) or separator.act(canonical(adj(raised)), vector):
            raise AssertionError('Spin partner failed to annihilate the maximum-spin multiplet')
        controls.append({'spin_z': str(F(3)-step), 'generated_amplitudes': len(vector), 'both_components_and_adjoints_annihilate': True})
        vector = separator.act(minus, vector)
    if vector:
        raise AssertionError('Maximum-spin ladder failed to terminate')
    singlet = {sum(1 << i for i in range(6)): F(1)}
    if separator.act(z, singlet) or separator.act(plus, singlet) or separator.act(minus, singlet):
        raise AssertionError('Closed-shell control is not a singlet')
    squares = sum(sum(a*a for a in separator.act(op, singlet).values())
                  for op in (C, canonical(adj(C)), raised, canonical(adj(raised))))
    direct = sum(singlet.get(state, F(0))*a for state, a in separator.act(P_scalar, singlet).items())
    if direct != squares or direct <= 0:
        raise AssertionError('Spin scalar failed its nonzero singlet positivity control')
    trace_moment = sum(c*seed(w, 12, 6) for w, c in P_scalar.items())
    threshold = -scalar_moment/(trace_moment-scalar_moment)
    return receipt | {'source_separator_sha256': cert['source_separator_sha256'],
        'spin_half_projection_exact': True, 'spin_scalar_completion_exact': True,
        'spin_scalar_polynomial_terms': len(P_scalar),
        'spin_scalar_moment': str(scalar_moment), 'spin_scalar_moment_float': float(scalar_moment),
        'spin_scalar_uniform_N6_trace_moment': str(trace_moment),
        'spin_scalar_strict_separation_below_trace_fraction': str(threshold),
        'spin_scalar_trace_threshold_float': float(threshold),
        'multiplet_action_controls': controls,
        'generated_amplitudes_in_controls': sum(row['generated_amplitudes'] for row in controls),
        'singlet_physical_control': {'generated_determinants': 1, 'exact_positive_expectation': str(direct),
            'independent_action_agrees_with_four_norm_squares': True},
        'full_fermionic_sector_basis_enumerated': False,
        'wall_seconds': time.monotonic()-start,
        'scope': 'A symmetry projection and spin-scalar completion of the established T1 law. Both reject the same frozen functional; neither is an independent molecular transfer or a new positivity law.'}


def run(replay=False):
    start = time.monotonic(); OUT.mkdir(exist_ok=True)
    gap_cert = json.loads((OUT/'gap_certificate.json').read_text())
    rows = [molecular_case('h6', 10, gap_cert), molecular_case('h8', 14)]
    source_bytes = (ROOT/'results/response_consistency_20260913/separator.json').read_bytes()
    parent_bytes = (ROOT/'results/trace_pricing_20260913/full_dual/witness.json').read_bytes()
    if replay:
        cert = json.loads((OUT/'spin_projected_separator.json').read_text())
    else:
        cert = projected_separator(source_bytes, parent_bytes)
        (OUT/'spin_projected_separator.json').write_text(json.dumps(cert, indent=2)+'\n')
        h6 = rows[0]; selected = next(row for row in h6['complement_gap_bounds'] if row['one_body_plus_R_nonnegative_already_suffices'])
        sector_cert = {'kind': 'molecular_double_occupancy_gap_v1', 'fixture_sha256': h6['fixture_sha256'],
            'tail_sha256': h6['tail_sha256'], 'doubly_occupied_spatial_orbital': selected['doubly_occupied_spatial_orbital'],
            'target_Ha': h6['complement_gap_target_Ha'],
            'delta_Ha': selected['shifted_D_lower_without_commutator_gap_Ha']}
        (OUT/'sector_gap_certificate.json').write_text(json.dumps(sector_cert, indent=2)+'\n')
    projected = check_projected(cert, source_bytes, parent_bytes)
    forbidden = [name for name in ('numpy', 'scipy', 'cvxpy', 'pyscf') if name in sys.modules]
    if forbidden:
        raise AssertionError('Numerical imports on the exact molecular diagnostic path')
    result = {'molecules': rows, 'spin_projected_separator': projected, 'numerical_packages_loaded': forbidden,
        'wall_seconds': time.monotonic()-start}
    (OUT/('molecular_replay.json' if replay else 'molecular_discovery.json')).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'molecules': len(rows), 'spin_scalar_moment': projected['spin_scalar_moment_float'],
        'wall_seconds': result['wall_seconds']}), flush=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(); parser.add_argument('--replay', action='store_true')
    run(parser.parse_args().replay)
