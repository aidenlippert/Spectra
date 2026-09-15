from pathlib import Path
from fractions import Fraction as F
import json
import pytest

from experiments.marginal_range_two_family_limit import replay


def certificate():
    return {'kind':'joint_diagonal_range2_family_limit_v4',
            'half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},
            'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4',
            'diagonal_shapes':[{102:1,612:-1}],'W':'1/10',
            'range_two_density_profile':['13/8','-11/8','-11/8','13/8'],
            'mixture':[{'weight':'1/3','vector':{0:7}},{'weight':'2/3','vector':{4095:11}}]}


def test_all_pair_moments_and_independent_vacuum_full_energy():
    result=replay(certificate())
    assert F(result['periodic_family_upper'])==F(13,5)
    assert result['charge_square_pair_moments']=={label:'0' for label in ['0,1','0,2','0,3','1,2']}
    assert all(v=='0' for v in result['quadratic_charge_moments'].values())


def test_prior_quadratic_dual_fails_new_consistency_constraint():
    root=Path(__file__).resolve().parents[1]
    c=json.loads((root/'results/marginal_graded_hubbard8/quadratic_charge_telescope/W_zero/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v4'
    with pytest.raises(ValueError,match='charge-square pair expectation'):replay(c)


def test_four_new_rows_bound_source_count_and_preserve_old_cap():
    c=certificate()
    c['mixture']=[{'weight':'1/16','vector':{0:i+1}} for i in range(16)]
    assert replay(c)['nearest_constraint_replay']['mixture_sources']==16
    c['kind']='joint_diagonal_range2_family_limit_v3'
    with pytest.raises(ValueError,match='scalar constraint'):replay(c)
    c['kind']='joint_diagonal_range2_family_limit_v4'
    c['mixture']=[{'weight':'1/17','vector':{0:i+1}} for i in range(17)]
    with pytest.raises(ValueError,match='scalar constraint'):replay(c)
