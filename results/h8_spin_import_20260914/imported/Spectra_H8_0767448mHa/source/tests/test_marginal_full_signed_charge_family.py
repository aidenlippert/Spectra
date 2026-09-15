from pathlib import Path
from fractions import Fraction as F
import json
import pytest
from experiments.marginal_range_two_family_limit import replay


def certificate():
    return {'kind':'joint_diagonal_range2_family_limit_v6','half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4','diagonal_shapes':[{102:1,612:-1}],'W':'1/10','range_two_density_profile':['13/8','-11/8','-11/8','13/8'],'mixture':[{'weight':'1/3','vector':{0:7}},{'weight':'2/3','vector':{4095:11}}]}


def test_all_signed_charge_moments_and_independent_physical_energy():
    result=replay(certificate())
    assert F(result['periodic_family_upper'])==F(13,5)
    assert len(result['signed_charge_moments'])==52
    assert all(v=='0' for v in result['signed_charge_moments'].values())


def test_prior_classically_extendible_indicator_dual_fails_signed_condition():
    root=Path(__file__).resolve().parents[1]
    c=json.loads((root/'results/marginal_graded_hubbard8/full_charge_indicators/W_zero/polished/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v6'
    with pytest.raises(ValueError,match='signed-charge pattern expectation'):replay(c)


def test_new_constraint_count_and_preserved_older_cap():
    c=certificate();c['mixture']=[{'weight':'1/58','vector':{0:i+1}} for i in range(58)]
    assert replay(c)['nearest_constraint_replay']['mixture_sources']==58
    c['kind']='joint_diagonal_range2_family_limit_v5'
    with pytest.raises(ValueError,match='scalar constraint'):replay(c)
    c['kind']='joint_diagonal_range2_family_limit_v6'
    c['mixture']=[{'weight':'1/59','vector':{0:i+1}} for i in range(59)]
    with pytest.raises(ValueError,match='scalar constraint'):replay(c)
