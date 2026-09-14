"""Exact fixed-profile range-two dual limits, reusing all nearest constraints."""
from fractions import Fraction as F

from experiments.marginal_diagonal_family_limit import replay as replay_nearest
from experiments.marginal_joint_family_limit import _physical_vector, _weight
from experiments.marginal_local_hubbard_block import _exact
from experiments.marginal_range_two_density import local_profile, diagonal_value


def replay(c):
    if type(c) is not dict or c.get('kind') != 'joint_diagonal_range2_family_limit_v1':
        raise ValueError('Explicit fixed-profile range-two family required')
    W = _exact(c.get('W'))
    profile = local_profile(6, W, c.get('range_two_density_profile'))
    nearest = dict(c, kind='joint_diagonal_family_limit_v1')
    nearest.pop('proposed_periodic_family_upper', None)
    nearest.pop('W', None)
    nearest.pop('range_two_density_profile', None)
    base = replay_nearest(nearest)
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
    return {'accepted': True, 'periodic_family_upper': str(upper),
            'nearest_constraint_replay': base, 'W': str(W),
            'range_two_density_profile': list(map(str, profile)),
            'range_two_local_expectation': str(extra),
            'scope': 'Upper limit on attainable periodic LOWER certificates for fixed U4,t1,V1/2,W, fixed range-two local density profile, projector sources/ratio/ceilings and supplied diagonal span. Six nearest-profile parameters and all supplied shape coefficients may vary. The range-two profile is fixed, not optimized or canceled by the dual. Not a physical ground-energy upper, general representability result or limit for other supports.'}
