"""Exact enlarged-family constraints and unchanged older refusal/cap gates."""
from pathlib import Path
import copy
import json
import pytest
from experiments.marginal_range_two_family_limit import replay
from experiments.marginal_diagonal_family_limit import _replay
from experiments.marginal_pair_transfer import LABELS
from test_marginal_spectator_hopping import family


def test_stationary_family_accepted_and_every_pair_moment_checked():
    c = family(); c['kind']='joint_diagonal_range2_family_limit_v10'
    r = replay(c)
    assert r['accepted'] and r['pair_transfer_moments'] == dict.fromkeys(LABELS,'0')
    assert r['periodic_family_upper'] == '13/5'
    assert len(r['spectator_hopping_moments']) == 14
    assert len(r['signed_moments'] if 'signed_moments' in r else r['signed_charge_moments']) == 52


@pytest.mark.parametrize('case', ['W_zero','W_plus_1'])
def test_preceding_violated_mixture_cannot_be_relabelled(case):
    root = Path(__file__).resolve().parents[1]
    c = json.loads((root/f'results/marginal_graded_hubbard8/spectator_hopping/{case}/thermal/final/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v10'
    with pytest.raises(ValueError,match='pair-transfer expectation'): replay(c)


def flags():
    return dict(free_range_two_profile=True,full_quadratic_charge=True,
                charge_square_pairs=True,full_charge_indicators=True,
                full_signed_charge=True,hopping_telescope=True,
                spin_telescope=True,spectator_hopping=True)


def test_new_source_cap_and_old_cap_are_distinct():
    c = family(); c.pop('W'); c.pop('range_two_density_profile')
    c['kind']='joint_diagonal_family_limit_v1'
    # One sparse shape gives77 old columns and81 with the four new constraints.
    c['mixture']=[{'weight':'1/81','vector':{0:1}}]*81
    with pytest.raises(ValueError,match='constraint count'): _replay(c,**flags())
    assert _replay(c,**flags(),pair_transfer=True)['mixture_sources'] == 81
    c['mixture']=[{'weight':'1/82','vector':{0:1}}]*82
    with pytest.raises(ValueError,match='constraint count'): _replay(c,**flags(),pair_transfer=True)


@pytest.mark.parametrize('pair_mode,spectator', [(1,True),(True,False)])
def test_explicit_hierarchy_required(pair_mode,spectator):
    c = family(); c.pop('W'); c.pop('range_two_density_profile')
    c['kind']='joint_diagonal_family_limit_v1'
    f=flags(); f['spectator_hopping']=spectator
    with pytest.raises(ValueError,match='requires spectator'): _replay(c,**f,pair_transfer=pair_mode)


def test_fixed_pair_coefficients_still_refused():
    c=family();c['kind']='joint_diagonal_range2_family_limit_v10';c['pair_transfer']={'0,1':'1/10'}
    with pytest.raises(ValueError,match='constrain moments'): replay(c)
