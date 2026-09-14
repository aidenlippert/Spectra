from pathlib import Path
from fractions import Fraction as F
import json
import pytest
from experiments.marginal_range_two_family_limit import replay


def certificate():
    return {'kind':'joint_diagonal_range2_family_limit_v5',
            'half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},
            'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4',
            'diagonal_shapes':[{102:1,612:-1}],'W':'1/10',
            'range_two_density_profile':['13/8','-11/8','-11/8','13/8'],
            'mixture':[{'weight':'1/3','vector':{0:7}},{'weight':'2/3','vector':{4095:11}}]}


def test_entire_indicator_basis_and_independent_physical_energy():
    result=replay(certificate())
    assert F(result['periodic_family_upper'])==F(13,5)
    assert len(result['charge_indicator_moments'])==12
    assert all(v=='0' for v in result['charge_indicator_moments'].values())


def test_previously_accepted_pair_dual_fails_higher_indicator_condition():
    root=Path(__file__).resolve().parents[1]
    c=json.loads((root/'results/marginal_graded_hubbard8/charge_square_pairs/W_zero/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v5'
    with pytest.raises(ValueError,match='charge indicator expectation'):replay(c)


def test_six_additional_rows_bound_sources_without_changing_old_cap():
    c=certificate();c['mixture']=[{'weight':'1/22','vector':{0:i+1}} for i in range(22)]
    assert replay(c)['nearest_constraint_replay']['mixture_sources']==22
    c['kind']='joint_diagonal_range2_family_limit_v4'
    with pytest.raises(ValueError,match='scalar constraint'):replay(c)
    c['kind']='joint_diagonal_range2_family_limit_v5'
    c['mixture']=[{'weight':'1/23','vector':{0:i+1}} for i in range(23)]
    with pytest.raises(ValueError,match='scalar constraint'):replay(c)
