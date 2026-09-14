from fractions import Fraction as F
from pathlib import Path
import json
import pytest
from experiments.marginal_spin_telescope import LABELS,coefficients,actions
from experiments.marginal_hopping_telescope import projected_matrix
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_charged_projectors import _particlehole,_spinflip
from test_marginal_hopping_telescope import toy


def direct(s,terms):
    image={};z=[((s>>(2*i))&1)-((s>>(2*i+1))&1) for i in range(6)]
    for label,a in terms.items():
        i,j=LABELS[label]
        for left,right,sign in [(i,j,1),(4-j,4-i,-1),(i+1,j+1,-1),(5-j,5-i,1)]:
            image[s]=image.get(s,0)+a*sign*z[left]*z[right]
            if z[left]*z[right]==-1:
                t=s^(3<<(2*left))^(3<<(2*right));image[t]=image.get(t,0)+2*a*sign
    return {t:a for t,a in image.items() if a}


def test_all_spin_components_against_independent_local_spin_actions():
    for terms in [{key:1} for key in LABELS]+[{key:F(i+1,7) for i,key in enumerate(LABELS)}]:
        action=actions(terms)
        for s in range(4096):
            assert action[s]==direct(s,terms)
            for transform in (_particlehole,_spinflip):
                us,phase=transform(s);expected={}
                for t,a in action[s].items():
                    ut,sign=transform(t);expected[ut]=sign*a
                assert expected=={t:phase*a for t,a in action[us].items()}


def test_projected_spin_image_has_no_discarded_components():
    terms={key:F(i+1,7) for i,key in enumerate(LABELS)};action=actions(terms)
    for _,_,cols in sector_matrices(6,0,0).values():
        matrix=projected_matrix(cols,action)
        for j,col in enumerate(cols):
            image={}
            for s,a in col.items():
                for t,b in direct(s,terms).items():image[t]=image.get(t,0)+a*b
            rebuilt={s:F(matrix[i][j]*a,sum(v*v for v in basis.values())) for i,basis in enumerate(cols) for s,a in basis.items()}
            assert {s:a for s,a in image.items() if a}=={s:a for s,a in rebuilt.items() if a}


def test_spin_pair_translates_cancel_exactly():
    for sites in (8,9,10,17):
        for i,j in LABELS.values():
            totals={}
            for offset in range(sites):
                for left,right,a in [(i,j,1),(4-j,4-i,-1),(i+1,j+1,-1),(5-j,5-i,1)]:
                    pair=tuple(sorted(((left+offset)%sites,(right+offset)%sites)));totals[pair]=totals.get(pair,0)+a
            assert not any(totals.values())


def test_exact_energy_accepts_safe_refuses_unsafe_and_old_version():
    c=toy();c.update(kind='hubbard_projector_extension_v12',spin_telescope={key:'1/100' for key in LABELS})
    r=replay(c)
    assert r['accepted'] and r['spin_telescope']==c['spin_telescope']
    assert r['target']==c['target'] and r['open_lower_density']==r['periodic_lower_density']=='-9/5'
    assert r['local_sum_dimensions']==4096 and r['local_maximum_psd_dimension']==200
    c['penalized_lower']='-7'
    with pytest.raises(ValueError):replay(c)
    c['kind']='hubbard_projector_extension_v11'
    with pytest.raises(ValueError,match='requires v12'):replay(c)


@pytest.mark.parametrize('source',[{}, {'0,0':1},{'3,4':1},{'0,1':0},{'0,1':True},{'0,1':.1},{'0,1':'1/1000001'},['0,1']])
def test_spin_coefficient_refusals(source):
    with pytest.raises(ValueError):coefficients(source)


def family():
    return {'kind':'joint_diagonal_range2_family_limit_v8','half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4','diagonal_shapes':[{102:1,612:-1}],'W':'1/10','range_two_density_profile':['1/8']*4,'mixture':[{'weight':'1/3','vector':{0:7}},{'weight':'2/3','vector':{4095:11}}]}


def test_spin_family_acceptance_and_previous_violated_witness_refusal():
    from experiments.marginal_range_two_family_limit import replay as family_replay
    c=family();r=family_replay(c)
    assert r['spin_telescope_moments']==dict.fromkeys(LABELS,'0')
    assert r['periodic_family_upper']=='13/5'
    old=json.loads((Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/hopping_telescope/W_zero/range_two_family_limit_certificate.json').read_text())
    old['kind']='joint_diagonal_range2_family_limit_v8'
    with pytest.raises(ValueError,match='spin telescope expectation'):family_replay(old)
    c['spin_telescope']={'0,1':'1/2'}
    with pytest.raises(ValueError,match='constrain moments'):family_replay(c)


def test_spin_family_mode_and_old_source_cap_preserved():
    from experiments.marginal_diagonal_family_limit import _replay
    c=family();c.pop('W');c.pop('range_two_density_profile');c['kind']='joint_diagonal_family_limit_v1';c['mixture']=[{'weight':'1/63','vector':{0:1}}]*63
    flags=dict(free_range_two_profile=True,full_quadratic_charge=True,charge_square_pairs=True,full_charge_indicators=True,full_signed_charge=True,hopping_telescope=True)
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**flags)
    assert _replay(c,**flags,spin_telescope=True)['mixture_sources']==63
    flags['hopping_telescope']=False
    with pytest.raises(ValueError,match='requires hopping'):_replay(c,**flags,spin_telescope=True)
    flags['hopping_telescope']=True
    with pytest.raises(ValueError,match='requires hopping'):_replay(c,**flags,spin_telescope=1)
