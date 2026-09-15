from fractions import Fraction as F
import pytest

from experiments.marginal_quadratic_charge_telescope import coefficients, five_site_value, local_value, PAIRS
from experiments.marginal_projector_extendibility import replay


def direct(s):
    q = [((s // 4**i) % 4 == 3)-((s // 4**i) % 4 == 0) for i in range(6)]
    return F(q[0]*q[3]-2*q[1]*q[4]+q[2]*q[5], 7)


def certificate():
    minimum = min(direct(s)+int(s % 1024 == 102)-int(s % 1024 == 612)
                  -int(s // 4 == 102)+int(s // 4 == 612) for s in range(4096))
    return {'kind': 'hubbard_projector_extension_v7', 'chain_sites': 10,
            'target': {'U': '0', 't': '0', 'V': '0', 'W': '0'},
            'local_window': {'kind': 'local_hubbard_range2_block_v1', 'sites': 6, 'U': '0', 't': '0', 'V': '0'},
            'vector': {0x999: 1, 0x666: 1}, 'windows': 2, 'projector_sum_ceiling': '2',
            'penalty': '0', 'penalized_lower': str(minimum),
            'joint': {'vector': {62: 1, 3008: 1}, 'windows': 2, 'ratio': '1/2', 'projector_sum_ceiling': '2', 'penalty': '0'},
            'telescoping_diagonal': {102: 1, 612: -1}, 'quadratic_charge_telescope': {'0,3': '1/7'}}


def test_local_full_fock_psd_matches_independent_three_pair_energy():
    c = certificate()
    result = replay(c)
    assert result['local_sum_dimensions'] == 4096
    assert result['local_maximum_psd_dimension'] == 200
    assert F(result['periodic_lower_density']) == F(c['penalized_lower'])/5
    assert result['quadratic_components'] == 1
    assert all(local_value(s, {'0,3': F(1, 7)}) == direct(s) for s in range(4096))
    c['penalized_lower'] = str(F(c['penalized_lower'])+F(1, 7))
    with pytest.raises(ValueError):
        replay(c)


def test_all_quadratic_components_are_reflection_odd_and_periodic_sum_is_zero():
    assert len(PAIRS) == 6
    terms = coefficients({f'{i},{j}': index+1 for index, (i, j) in enumerate(PAIRS)})
    for s in range(1024):
        reflected = sum(((s >> (2*i)) & 3) << (2*(4-i)) for i in range(5))
        assert five_site_value(reflected, terms) == -five_site_value(s, terms)
    for s in range(4**8):
        q = [((s >> (2*i)) & 3).bit_count()-1 for i in range(8)]
        assert sum(q[i]*q[(i+3) % 8]-2*q[(i+1) % 8]*q[(i+4) % 8]+q[(i+2) % 8]*q[(i+5) % 8] for i in range(8)) == 0


@pytest.mark.parametrize('terms', [{}, {'0,4': 1}, {'3,0': 1}, {'1,4': 1}, {'0,3': 0}, {'0,3': 0.1}, {'0,3': '1/1000001'}, {'0,3': True}])
def test_malformed_or_noncanonical_compact_terms_refused(terms):
    with pytest.raises(ValueError):
        coefficients(terms)


def test_old_energy_version_refuses_compact_quadratic_field():
    c = certificate()
    c['kind'] = 'hubbard_projector_extension_v6'
    with pytest.raises(ValueError, match='requires v7'):
        replay(c)
