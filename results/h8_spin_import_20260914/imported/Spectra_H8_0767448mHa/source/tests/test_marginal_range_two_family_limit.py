from fractions import Fraction as F
import pytest

from experiments.marginal_range_two_family_limit import replay


def certificate(W):
    return {'kind': 'joint_diagonal_range2_family_limit_v1',
            'half_vector': {0x999: 1, 0x666: 1}, 'charged_vector': {62: 1, 3008: 1},
            'ratio': '1/2', 'theta_half': '1/2', 'theta_joint': '3/4',
            'diagonal_shapes': [{102: 1, 612: -1}], 'W': str(W),
            'range_two_density_profile': [str(F(5, 4)*W)]*4,
            'mixture': [{'weight': '1/3', 'vector': {0: 7}}, {'weight': '2/3', 'vector': {4095: 11}}]}


@pytest.mark.parametrize('W', [F(-1, 5), F(0), F(1, 10)])
def test_independent_vacuum_full_energy_for_fixed_range_two_profile(W):
    result = replay(certificate(W))
    assert F(result['periodic_family_upper']) == F(5, 2)+W
    assert F(result['range_two_local_expectation']) == 5*W
    assert result['nearest_constraint_replay']['diagonal_moments'] == ['0']


def test_wrong_profile_claim_and_inherited_constraint_refusals():
    for field, value in [('range_two_density_profile', ['0']*4),
                         ('proposed_periodic_family_upper', '0'), ('W', 0.1),
                         ('kind', 'joint_diagonal_family_limit_v1'),
                         ('mixture', [{'weight': '1/2', 'vector': {0: 1}}])]:
        c = certificate(F(1, 10))
        c[field] = value
        with pytest.raises(ValueError):
            replay(c)


def test_nearest_only_verifier_refuses_range_two_fields():
    from experiments.marginal_diagonal_family_limit import replay as replay_nearest
    c = certificate(F(1, 10))
    c['kind'] = 'joint_diagonal_family_limit_v1'
    with pytest.raises(ValueError, match='Range-two fields'):
        replay_nearest(c)
