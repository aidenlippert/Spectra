from fractions import Fraction as F
from pathlib import Path
import json
import pytest

from experiments.marginal_range_two_family_limit import replay


def certificate():
    return {'kind': 'joint_diagonal_range2_family_limit_v2',
            'half_vector': {0x999: 1, 0x666: 1}, 'charged_vector': {62: 1, 3008: 1},
            'ratio': '1/2', 'theta_half': '1/2', 'theta_joint': '3/4',
            'diagonal_shapes': [{102: 1, 612: -1}], 'W': '1/10',
            'range_two_density_profile': ['1/8']*4,
            'mixture': [{'weight': '1/3', 'vector': {0: 7}}, {'weight': '2/3', 'vector': {4095: 11}}]}


def test_free_profile_cap_is_invariant_under_large_signed_profile_changes():
    c = certificate()
    uniform = replay(c)
    c['range_two_density_profile'] = ['13/8', '-11/8', '-11/8', '13/8']
    changed = replay(c)
    assert uniform['periodic_family_upper'] == changed['periodic_family_upper'] == '13/5'
    assert changed['range_two_profile_gradient'] == '0'


def test_previously_accepted_fixed_profile_dual_fails_the_new_constraint():
    root = Path(__file__).resolve().parents[1]
    c = json.loads((root/'results/marginal_graded_hubbard8/range_two_density/W_plus_1_10/range_two_family_limit_certificate.json').read_text())
    c['kind'] = 'joint_diagonal_range2_family_limit_v2'
    with pytest.raises(ValueError, match='profile expectation'):
        replay(c)


def test_one_additional_constraint_allows_only_one_additional_source():
    c = certificate()
    c['mixture'] = [{'weight': '1/11', 'vector': {0: i+1}} for i in range(11)]
    assert replay(c)['nearest_constraint_replay']['mixture_sources'] == 11
    c['kind'] = 'joint_diagonal_range2_family_limit_v1'
    with pytest.raises(ValueError, match='scalar constraint'):
        replay(c)
    c['kind'] = 'joint_diagonal_range2_family_limit_v2'
    c['mixture'] = [{'weight': '1/12', 'vector': {0: i+1}} for i in range(12)]
    with pytest.raises(ValueError, match='scalar constraint'):
        replay(c)
