from fractions import Fraction as F
import copy
import pytest

from experiments.marginal_range_two_density import local_profile, diagonal_value
from experiments.marginal_projector_extendibility import replay


def certificate(profile):
    y = {102: F(1), 612: F(-1)}
    def value(s):
        digits = [(s // 4**i) % 4 for i in range(6)]
        q = [(d == 3)-(d == 0) for d in digits]
        return sum(profile[i]*q[i]*q[i+2] for i in range(4))+y.get(s % 1024, 0)-y.get(s // 4, 0)
    lower = min(value(s) for s in range(4096))
    return {'kind': 'hubbard_projector_extension_v6', 'chain_sites': 10,
            'target': {'U': '0', 't': '0', 'V': '0', 'W': '1/5'},
            'local_window': {'kind': 'local_hubbard_range2_block_v1', 'sites': 6,
                             'U': '0', 't': '0', 'V': '0', 'range_two_density_profile': list(map(str, profile))},
            'vector': {0x999: 1, 0x666: 1}, 'windows': 2, 'projector_sum_ceiling': '2',
            'penalty': '0', 'penalized_lower': str(lower),
            'joint': {'vector': {62: 1, 3008: 1}, 'windows': 2, 'ratio': '1/2',
                      'projector_sum_ceiling': '2', 'penalty': '0'},
            'telescoping_diagonal': y}


@pytest.mark.parametrize('profile', [[F(1, 10), F(2, 5), F(2, 5), F(1, 10)],
                                    [F(3, 5), F(-1, 10), F(-1, 10), F(3, 5)]])
def test_full_fock_lower_matches_independent_diagonal_minimum(profile):
    c = certificate(profile)
    result = replay(c)
    assert result['local_sum_dimensions'] == 4096
    assert result['local_maximum_psd_dimension'] == 200
    assert result['target']['W'] == '1/5'
    assert result['closing_interaction_norm'] == '2/5'
    assert F(result['periodic_lower_density']) == F(c['penalized_lower'])/5
    assert F(result['open_lower_density']) == F(c['penalized_lower'])/5-F(1, 25)
    c['penalized_lower'] = str(F(c['penalized_lower'])+F(1, 10))
    with pytest.raises(ValueError):
        replay(c)


def test_periodic_range_two_cover_and_two_opening_bonds_on_all_determinants():
    profile = local_profile(6, F(1, 5), ['1/10', '2/5', '2/5', '1/10'])
    integer_profile = [int(10*v) for v in profile]
    for s in range(4**8):
        digits = [(s // 4**i) % 4 for i in range(8)]
        q = [(d == 3)-(d == 0) for d in digits]
        translated = sum(integer_profile[j]*q[(i+j) % 8]*q[(i+j+2) % 8]
                         for i in range(8) for j in range(4))
        periodic = sum(q[i]*q[(i+2) % 8] for i in range(8))
        opened = sum(q[i]*q[i+2] for i in range(6))
        assert translated == 10*periodic
        assert abs(periodic-opened) <= 2


def test_nonzero_range_two_cannot_be_smuggled_into_old_version():
    c = certificate([F(1, 4)]*4)
    c['kind'] = 'hubbard_projector_extension_v5'
    c['local_window']['kind'] = 'local_hubbard_block_v1'
    with pytest.raises(ValueError, match='requires v6'):
        replay(c)


@pytest.mark.parametrize('field,value', [('W', None), ('W', 0.2),
                                        ('profile', ['1/5']*4), ('profile', ['1/10', '2/5', '1/10', '2/5']),
                                        ('profile', [1, 2]), ('profile', [0.25]*4),
                                        ('kind', 'local_hubbard_block_v1')])
def test_malformed_range_two_target_or_profile_refused(field, value):
    c = certificate([F(1, 4)]*4)
    if field == 'W':
        c['target']['W'] = value
    elif field == 'profile':
        c['local_window']['range_two_density_profile'] = value
    else:
        c['local_window']['kind'] = value
    with pytest.raises(ValueError):
        replay(c)
