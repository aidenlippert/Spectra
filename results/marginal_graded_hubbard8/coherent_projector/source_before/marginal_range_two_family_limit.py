"""Exact fixed-profile range-two dual limits, reusing all nearest constraints."""
from fractions import Fraction as F

from experiments.marginal_diagonal_family_limit import _replay as replay_constraints
from experiments.marginal_joint_family_limit import _physical_vector, _weight
from experiments.marginal_local_hubbard_block import _exact
from experiments.marginal_range_two_density import local_profile, diagonal_value


def replay(c):
    if type(c) is not dict or c.get('kind') not in ('joint_diagonal_range2_family_limit_v1', 'joint_diagonal_range2_family_limit_v2', 'joint_diagonal_range2_family_limit_v3', 'joint_diagonal_range2_family_limit_v4', 'joint_diagonal_range2_family_limit_v5', 'joint_diagonal_range2_family_limit_v6', 'joint_diagonal_range2_family_limit_v7', 'joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12'):
        raise ValueError('Explicit range-two family version required')
    free_profile = c['kind'] != 'joint_diagonal_range2_family_limit_v1'
    full_quadratic = c['kind'] in ('joint_diagonal_range2_family_limit_v3', 'joint_diagonal_range2_family_limit_v4', 'joint_diagonal_range2_family_limit_v5', 'joint_diagonal_range2_family_limit_v6', 'joint_diagonal_range2_family_limit_v7', 'joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12')
    square_pairs = c['kind'] in ('joint_diagonal_range2_family_limit_v4', 'joint_diagonal_range2_family_limit_v5', 'joint_diagonal_range2_family_limit_v6', 'joint_diagonal_range2_family_limit_v7', 'joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12')
    full_indicators = c['kind'] in ('joint_diagonal_range2_family_limit_v5', 'joint_diagonal_range2_family_limit_v6', 'joint_diagonal_range2_family_limit_v7', 'joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12')
    full_signed = c['kind'] in ('joint_diagonal_range2_family_limit_v6', 'joint_diagonal_range2_family_limit_v7', 'joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12')
    coherent = c['kind'] in ('joint_diagonal_range2_family_limit_v7', 'joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12')
    spin = c['kind'] in ('joint_diagonal_range2_family_limit_v8', 'joint_diagonal_range2_family_limit_v9', 'joint_diagonal_range2_family_limit_v10', 'joint_diagonal_range2_family_limit_v11', 'joint_diagonal_range2_family_limit_v12')
    three_spectator = c['kind'] == 'joint_diagonal_range2_family_limit_v12'
    two_spectator = three_spectator or c['kind'] == 'joint_diagonal_range2_family_limit_v11'
    pair_transfer = two_spectator or c['kind'] == 'joint_diagonal_range2_family_limit_v10'
    spectator = pair_transfer or c['kind'] == 'joint_diagonal_range2_family_limit_v9'
    W = _exact(c.get('W'))
    profile = local_profile(6, W, c.get('range_two_density_profile'))
    nearest = dict(c, kind='joint_diagonal_family_limit_v1')
    nearest.pop('proposed_periodic_family_upper', None)
    nearest.pop('W', None)
    nearest.pop('range_two_density_profile', None)
    base = replay_constraints(nearest, free_range_two_profile=free_profile, full_quadratic_charge=full_quadratic, charge_square_pairs=square_pairs, full_charge_indicators=full_indicators, full_signed_charge=full_signed, hopping_telescope=coherent, spin_telescope=spin, spectator_hopping=spectator, pair_transfer=pair_transfer, two_spectator_hopping=two_spectator, three_spectator_hopping=three_spectator)
    # Existing replay already proves trace, PSD, all six profile and all
    # selected telescoping moments, and both fixed fidelity inequalities.
    extra = F(0)
    for item in c['mixture']:
        vector, norm = _physical_vector(item['vector'])
        extra += _weight(item['weight'])*sum((a*a*diagonal_value(s, profile)
                                             for s, a in vector.items()), F(0))/norm
    upper = F(base['periodic_family_upper'])+extra/5
    if 'proposed_periodic_family_upper' in c and F(c['proposed_periodic_family_upper']) != upper:
        raise ValueError('Range-two family upper disagrees with physical expectation')
    result = {'accepted': True, 'periodic_family_upper': str(upper),
            'nearest_constraint_replay': base, 'W': str(W),
            'range_two_density_profile': list(map(str, profile)),
            'range_two_local_expectation': str(extra),
            'scope': 'Upper limit on attainable periodic LOWER certificates for fixed U4,t1,V1/2,W, fixed range-two local density profile, projector sources/ratio/ceilings and supplied diagonal span. Six nearest-profile parameters and all supplied shape coefficients may vary. The range-two profile is fixed, not optimized or canceled by the dual. Not a physical ground-energy upper, general representability result or limit for other supports.'}
    if free_profile:
        result['range_two_profile_gradient'] = base['range_two_profile_gradient']
        result['scope'] = 'Upper limit on attainable periodic LOWER certificates for fixed U4,t1,V1/2,W, projector sources/ratio/ceilings and supplied diagonal span. All reflected mean-correct nearest and range-two density profiles are covered: the additional exact [1,-1,-1,1] moment is zero. All supplied shape coefficients may vary. Not a physical ground-energy upper or a cap over other supports or sources.'
    if full_quadratic:
        result['quadratic_charge_moments'] = base['quadratic_charge_moments']
        result['scope'] += ' All six reflection-odd quadratic five-site charge telescopes are also covered by exact zero moments.'
    if square_pairs:
        result['charge_square_pair_moments'] = base['charge_square_pair_moments']
        result['scope'] += ' All four reflection-odd pairs of charge squares are covered by exact zero moments.'
    if full_indicators:
        result['charge_indicator_moments'] = base['charge_indicator_moments']
        result['scope'] += ' The entire twelve-dimensional reflection-odd space of functions of five binary empty/double indicators is covered by exact zero moments.'
    if full_signed:
        result['signed_charge_moments'] = base['signed_charge_moments']
        result['scope'] += ' All 52 PH-even reflection-odd functions of five signed charges are covered by exact zero moments.'
    if coherent:
        result['hopping_telescope_moment'] = base['hopping_telescope_moment']
        result['scope'] += ' The unrestricted range-three hopping telescope coefficient is covered by its exact zero moment.'
    if spin:
        result['spin_telescope_moments'] = base['spin_telescope_moments']
        result['scope'] += ' All four reflection-odd spin-dot telescopes are covered by exact zero moments.'
    if spectator:
        result['spectator_hopping_moments'] = base['spectator_hopping_moments']
        result['scope'] += ' All fourteen one-spectator charge-conditioned hopping telescopes are covered by exact zero moments.'
    if pair_transfer:
        result['pair_transfer_moments'] = base['pair_transfer_moments']
        result['scope'] += ' All four pair-transfer telescopes are covered by exact zero moments.'
    if two_spectator:
        result['two_spectator_hopping_moments'] = base['two_spectator_hopping_moments']
        result['scope'] += ' All thirty two-spectator charge-hopping telescopes are covered by exact zero moments.'
    if three_spectator:
        result['three_spectator_hopping_moments'] = base['three_spectator_hopping_moments']
    return result
