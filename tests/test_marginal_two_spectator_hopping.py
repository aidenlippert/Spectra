"""Exact two-spectator energy/family extension and preserved refusal gates."""
from fractions import Fraction as F
from pathlib import Path
import json
import pytest
from experiments.marginal_two_spectator_hopping import LABELS,coefficients,actions
from experiments.marginal_hopping_telescope import projected_matrix
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_range_two_family_limit import replay as family_replay
from experiments.marginal_diagonal_family_limit import _replay
from test_marginal_two_spectator_probe import module as independent
from test_marginal_pair_transfer import pair_toy
from test_marginal_spectator_hopping import family
from test_marginal_pair_family import flags


def test_every_production_action_matches_the_independent_full_fock_probe():
    assert set(LABELS.values())==set(independent.LABELS)
    for key,label in LABELS.items():assert actions({key:1})==independent.operator(label)


def test_mixed_projection_reconstructs_entire_physical_image():
    terms={key:F(i+1,101) for i,key in enumerate(LABELS)};data=actions(terms)
    for _,_,cols in sector_matrices(6,0,0).values():
        matrix=projected_matrix(cols,data)
        for j,col in enumerate(cols):
            image={}
            for s,a in col.items():
                for t,b in data[s].items():image[t]=image.get(t,0)+a*b
            rebuilt={s:F(matrix[i][j]*a,sum(v*v for v in basis.values()))
                     for i,basis in enumerate(cols) for s,a in basis.items()}
            assert {s:a for s,a in image.items() if a}=={s:a for s,a in rebuilt.items() if a}


def two_toy():
    c=pair_toy();c['kind']='hubbard_projector_extension_v15'
    c['two_spectator_hopping']={key:'1/10000' for key in LABELS}
    return c


def test_safe_energy_and_negative_psd_refusal():
    c=two_toy();r=replay(c)
    assert r['accepted'] and r['two_spectator_hopping']==c['two_spectator_hopping']
    assert r['local_sum_dimensions']==4096 and len(r['local_sectors'])==94
    assert r['open_lower_density']==r['periodic_lower_density']=='-9/5'
    c['two_spectator_hopping']={next(iter(LABELS)):'100'}
    with pytest.raises(ValueError):replay(c)


@pytest.mark.parametrize('version',[1,2,4,5,6,7,8,9,10,11,12,13,14])
def test_old_versions_cannot_ignore_two_spectator_terms(version):
    c=two_toy();c['kind']=f'hubbard_projector_extension_v{version}'
    with pytest.raises(ValueError,match='requires v15'):replay(c)


@pytest.mark.parametrize('value',[None,{}, {'0,1,2,3,1,2':1},{'0,1,0,3,1,1':1},
    {'0,1,2,3,1,1':0},{'0,1,2,3,1,1':True},{'0,1,2,3,1,1':.1},
    {'0,1,2,3,1,1':'1/1000001'}])
def test_invalid_exact_coefficients_refused(value):
    with pytest.raises(ValueError):coefficients(value)


def test_stationary_family_has_all_thirty_zero_moments():
    c=family();c['kind']='joint_diagonal_range2_family_limit_v11';r=family_replay(c)
    assert r['two_spectator_hopping_moments']==dict.fromkeys(LABELS,'0')
    assert r['periodic_family_upper']=='13/5'
    assert len(r['pair_transfer_moments'])==4
    c['two_spectator_hopping']={next(iter(LABELS)):'1/1000'}
    with pytest.raises(ValueError,match='constrain moments'):family_replay(c)


@pytest.mark.parametrize('case',['W_zero','W_plus_1'])
def test_previous_violated_mixtures_cannot_be_relabelled(case):
    root=Path(__file__).resolve().parents[1]
    c=json.loads((root/f'results/marginal_graded_hubbard8/pair_transfer/{case}/family/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v11'
    with pytest.raises(ValueError,match='two-spectator hopping expectation'):family_replay(c)


def test_enlarged_source_limit_and_required_hierarchy():
    c=family();c.pop('W');c.pop('range_two_density_profile');c['kind']='joint_diagonal_family_limit_v1'
    f=flags();f['pair_transfer']=True
    c['mixture']=[{'weight':'1/111','vector':{0:1}}]*111
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**f)
    assert _replay(c,**f,two_spectator_hopping=True)['mixture_sources']==111
    c['mixture']=[{'weight':'1/112','vector':{0:1}}]*112
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**f,two_spectator_hopping=True)
    f['pair_transfer']=False
    with pytest.raises(ValueError,match='requires pair-transfer'):_replay(c,**f,two_spectator_hopping=True)
    f['pair_transfer']=True
    with pytest.raises(ValueError,match='requires pair-transfer'):_replay(c,**f,two_spectator_hopping=1)
