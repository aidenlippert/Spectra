from fractions import Fraction as F
from pathlib import Path
import json
import pytest
from experiments.marginal_spectator_hopping import LABELS,coefficients,actions
from experiments.marginal_hopping_telescope import projected_matrix
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_charged_projectors import _particlehole,_spinflip
from test_marginal_hopping_telescope import toy


def direct(state,terms):
    image={};q=[((state>>(2*i))&3).bit_count()-1 for i in range(6)]
    for label,value in terms.items():
        i,j,k,power=LABELS[label]
        for left,right,spectator,sign in [(i,j,k,1),(4-j,4-i,4-k,-1),(i+1,j+1,k+1,-1),(5-j,5-i,5-k,1)]:
            for spin in (0,1):
                a,b=2*left+spin,2*right+spin
                if ((state>>a)&1)==((state>>b)&1):continue
                phase=(-1)**sum((state>>m)&1 for m in range(a+1,b));target=state^(1<<a)^(1<<b)
                image[target]=image.get(target,0)+value*sign*q[spectator]**power*phase
    return {s:a for s,a in image.items() if a}


def test_complete_parity_selected_basis_and_independent_all_fock_actions():
    assert len(LABELS)==14 and sum(v[-1]==2 for v in LABELS.values())==9
    for terms in [{key:1} for key in LABELS]+[{key:F(i+1,101) for i,key in enumerate(LABELS)}]:
        data=actions(terms)
        for s in range(4096):
            assert data[s]==direct(s,terms)
            for transform in (_particlehole,_spinflip):
                us,phase=transform(s);right={}
                for t,a in data[s].items():ut,sign=transform(t);right[ut]=sign*a
                assert right=={t:phase*a for t,a in data[us].items()}


def test_spectator_projection_retains_the_entire_physical_image():
    terms={key:F(i+1,101) for i,key in enumerate(LABELS)};data=actions(terms)
    for _,_,cols in sector_matrices(6,0,0).values():
        matrix=projected_matrix(cols,data)
        for j,col in enumerate(cols):
            image={}
            for s,a in col.items():
                for t,b in direct(s,terms).items():image[t]=image.get(t,0)+a*b
            rebuilt={s:F(matrix[i][j]*a,sum(v*v for v in basis.values())) for i,basis in enumerate(cols) for s,a in basis.items()}
            assert {s:a for s,a in image.items() if a}=={s:a for s,a in rebuilt.items() if a}


def test_all_spectator_translates_cancel():
    for sites in (8,9,10,17):
        for i,j,k,power in LABELS.values():
            totals={}
            for offset in range(sites):
                for a,b,c,sign in [(i,j,k,1),(4-j,4-i,4-k,-1),(i+1,j+1,k+1,-1),(5-j,5-i,5-k,1)]:
                    pair=tuple(sorted(((a+offset)%sites,(b+offset)%sites)));key=pair+((c+offset)%sites,power);totals[key]=totals.get(key,0)+sign
            assert not any(totals.values())


def test_exact_energy_acceptance_negative_psd_and_old_version_refusal():
    c=toy();c.update(kind='hubbard_projector_extension_v13',spin_telescope={key:'1/100' for key in ['0,1','0,2','0,3','1,2']},spectator_hopping={key:'1/1000' for key in LABELS})
    r=replay(c)
    assert r['accepted'] and r['spectator_hopping']==c['spectator_hopping']
    assert r['target']==c['target'] and r['open_lower_density']==r['periodic_lower_density']=='-9/5'
    assert r['local_sum_dimensions']==4096 and r['local_maximum_psd_dimension']==200
    c['penalized_lower']='-7'
    with pytest.raises(ValueError):replay(c)
    c['kind']='hubbard_projector_extension_v12'
    with pytest.raises(ValueError,match='requires v13'):replay(c)


@pytest.mark.parametrize('source',[{}, {'0,1,0':1},{'3,4,2':1},{'0,1,2,2':1},{'0,1,2':0},{'0,1,2':True},{'0,1,2':.1},{'0,1,2':'1/1000001'}])
def test_bad_spectator_coefficients_refused(source):
    with pytest.raises(ValueError):coefficients(source)


def family():
    return {'kind':'joint_diagonal_range2_family_limit_v9','half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4','diagonal_shapes':[{102:1,612:-1}],'W':'1/10','range_two_density_profile':['1/8']*4,'mixture':[{'weight':'1/3','vector':{0:7}},{'weight':'2/3','vector':{4095:11}}]}


def test_spectator_family_acceptance_and_old_violated_witness_refusal():
    from experiments.marginal_range_two_family_limit import replay as family_replay
    c=family();r=family_replay(c)
    assert r['spectator_hopping_moments']==dict.fromkeys(LABELS,'0')
    assert r['periodic_family_upper']=='13/5'
    old=json.loads((Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/spin_telescope/W_zero/polished/range_two_family_limit_certificate.json').read_text())
    old['kind']='joint_diagonal_range2_family_limit_v9'
    with pytest.raises(ValueError,match='spectator hopping expectation'):family_replay(old)
    c['spectator_hopping']={'0,1,2':'1/2'}
    with pytest.raises(ValueError,match='constrain moments'):family_replay(c)


def test_spectator_family_mode_and_old_source_cap_preserved():
    from experiments.marginal_diagonal_family_limit import _replay
    c=family();c.pop('W');c.pop('range_two_density_profile');c['kind']='joint_diagonal_family_limit_v1';c['mixture']=[{'weight':'1/77','vector':{0:1}}]*77
    flags=dict(free_range_two_profile=True,full_quadratic_charge=True,charge_square_pairs=True,full_charge_indicators=True,full_signed_charge=True,hopping_telescope=True,spin_telescope=True)
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**flags)
    assert _replay(c,**flags,spectator_hopping=True)['mixture_sources']==77
    flags['spin_telescope']=False
    with pytest.raises(ValueError,match='requires spin'):_replay(c,**flags,spectator_hopping=True)
    flags['spin_telescope']=True
    with pytest.raises(ValueError,match='requires spin'):_replay(c,**flags,spectator_hopping=1)
