from fractions import Fraction as F
from pathlib import Path
import json
import pytest

from experiments.marginal_range_two_family_limit import replay


def certificate():
    return {'kind': 'joint_diagonal_range2_family_limit_v3',
            'half_vector': {0x999: 1, 0x666: 1}, 'charged_vector': {62: 1, 3008: 1},
            'ratio': '1/2', 'theta_half': '1/2', 'theta_joint': '3/4',
            'diagonal_shapes': [{102: 1, 612: -1}], 'W': '1/10',
            'range_two_density_profile': ['13/8', '-11/8', '-11/8', '13/8'],
            'mixture': [{'weight': '1/3', 'vector': {0: 7}}, {'weight': '2/3', 'vector': {4095: 11}}]}


def test_full_quadratic_family_checks_all_six_moments_and_physical_energy():
    result = replay(certificate())
    assert F(result['periodic_family_upper']) == F(13, 5)
    assert result['quadratic_charge_moments'] == {label: '0' for label in ['0,0', '0,1', '0,2', '0,3', '1,1', '1,2']}
    assert result['range_two_profile_gradient'] == '0'


def test_previously_accepted_free_profile_dual_fails_quadratic_constraint():
    root = Path(__file__).resolve().parents[1]
    c = json.loads((root/'results/marginal_graded_hubbard8/free_range_two_profile/W_zero/range_two_family_limit_certificate.json').read_text())
    c['kind'] = 'joint_diagonal_range2_family_limit_v3'
    with pytest.raises(ValueError, match='quadratic charge expectation'):
        replay(c)


def test_one_extra_independent_row_permits_only_one_extra_source():
    c = certificate()
    c['mixture'] = [{'weight': '1/12', 'vector': {0: i+1}} for i in range(12)]
    assert replay(c)['nearest_constraint_replay']['mixture_sources'] == 12
    c['kind'] = 'joint_diagonal_range2_family_limit_v2'
    with pytest.raises(ValueError, match='scalar constraint'):
        replay(c)
    c['kind'] = 'joint_diagonal_range2_family_limit_v3'
    c['mixture'] = [{'weight': '1/13', 'vector': {0: i+1}} for i in range(13)]
    with pytest.raises(ValueError, match='scalar constraint'):
        replay(c)
