from fractions import Fraction as F
import pytest
from experiments.marginal_hopping_telescope import actions,coefficient,projected_matrix
from experiments.marginal_local_hubbard_block import sector_matrices
from experiments.marginal_projector_extendibility import replay


def independent_action(s):
    out={}
    for i,weight in enumerate((1,-2,1)):
        for spin in (0,1):
            a,b=2*i+spin,2*(i+3)+spin
            if ((s>>a)&1)==((s>>b)&1):continue
            phase=(-1)**sum((s>>m)&1 for m in range(a+1,b))
            out[s^(1<<a)^(1<<b)]=weight*phase
    return out


def test_all_fock_car_and_projected_matrix_agree_with_bit_swap():
    action=actions()
    assert action==[independent_action(s) for s in range(4096)]
    data=sector_matrices(6,0,0)
    assert sum(len(cols) for _,_,cols in data.values())==4096
    for _,_,cols in data.values():
        matrix=projected_matrix(cols,action)
        lookup={s:(i,a) for i,col in enumerate(cols) for s,a in col.items()}
        for j,col in enumerate(cols):
            image={}
            for s,a in col.items():
                for t,b in independent_action(s).items():image[t]=image.get(t,0)+a*b
            # Reconstruct T E_j from the projected coordinates: verifies no
            # discarded components leak outside the supplied reflection block.
            rebuilt={s:F(matrix[i][j]*a,sum(v*v for v in basis.values())) for i,basis in enumerate(cols) for s,a in basis.items()}
            assert {s:a for s,a in image.items() if a}=={s:a for s,a in rebuilt.items() if a}


def test_periodic_car_coefficients_cancel():
    for sites in (8,9,10,17):
        totals={}
        for offset in range(sites):
            for i,weight in enumerate((1,-2,1)):
                for spin in (0,1):
                    a,b=2*((offset+i)%sites)+spin,2*((offset+i+3)%sites)+spin
                    for word in ((a,b),(b,a)):totals[word]=totals.get(word,0)+weight
        assert not any(totals.values())


def toy():
    return {'kind':'hubbard_projector_extension_v11','chain_sites':10,'target':{'U':'0','t':'0','V':'0','W':'0'},'local_window':{'kind':'local_hubbard_range2_block_v1','sites':6,'U':'0','t':'0','V':'0'},'vector':{0x999:1,0x666:1},'windows':2,'projector_sum_ceiling':'2','penalty':'0','penalized_lower':'-9','joint':{'vector':{62:1,3008:1},'windows':2,'ratio':'1/2','projector_sum_ceiling':'2','penalty':'0'},'telescoping_diagonal':{102:'1/1000',612:'-1/1000'},'quadratic_charge_telescope':{'0,3':'1/1000'},'charge_square_pair_telescope':{'0,1':'1/1000'},'higher_charge_indicator_telescope':{'0,1,2':'1/1000'},'signed_charge_telescope':{'-1,-1,-1,-1,0':'1/1000'},'hopping_telescope':'1'}


def test_exact_energy_accepts_safe_bound_refuses_unsafe_and_old_version():
    c=toy();r=replay(c)
    assert r['accepted'] and r['hopping_telescope']=='1'
    assert r['target']==c['target']
    assert r['open_lower_density']==r['periodic_lower_density']=='-9/5'
    assert r['local_sum_dimensions']==4096 and len(r['local_sectors'])==94
    assert r['local_maximum_psd_dimension']==200
    c['penalized_lower']='-7'
    with pytest.raises(ValueError):replay(c)
    c['kind']='hubbard_projector_extension_v10'
    with pytest.raises(ValueError,match='requires v11'):replay(c)


@pytest.mark.parametrize('value',[None,0,True,.1,'1/1000001','1000001',{},[1,-2,1]])
def test_hopping_coefficient_refusals(value):
    with pytest.raises(ValueError):coefficient(value)


def test_hopping_dual_accepts_stationary_mixture_refuses_previous_separator():
    import json
    from pathlib import Path
    from experiments.marginal_range_two_family_limit import replay as family_replay
    c={'kind':'joint_diagonal_range2_family_limit_v7','half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4','diagonal_shapes':[{102:1,612:-1}],'W':'1/10','range_two_density_profile':['1/8']*4,'mixture':[{'weight':'1/3','vector':{0:7}},{'weight':'2/3','vector':{4095:11}}]}
    result=family_replay(c)
    assert result['hopping_telescope_moment']=='0'
    assert result['periodic_family_upper']=='13/5'
    assert len(result['signed_charge_moments'])==52
    root=Path(__file__).resolve().parents[1]
    old=json.loads((root/'results/marginal_graded_hubbard8/full_signed_charge/W_zero/range_two_family_limit_certificate.json').read_text())
    old['kind']='joint_diagonal_range2_family_limit_v7'
    with pytest.raises(ValueError,match='Hopping telescope expectation'):family_replay(old)
    c['hopping_telescope']='1/2'
    with pytest.raises(ValueError,match='constrain moments'):family_replay(c)


def test_hopping_dual_requires_full_mode_and_preserves_old_source_cap():
    from experiments.marginal_diagonal_family_limit import _replay
    c={'kind':'joint_diagonal_family_limit_v1','half_vector':{0x999:1,0x666:1},'charged_vector':{62:1,3008:1},'ratio':'1/2','theta_half':'1/2','theta_joint':'3/4','diagonal_shapes':[{102:1,612:-1}],'mixture':[{'weight':'1/59','vector':{0:1}}]*59}
    flags=dict(free_range_two_profile=True,full_quadratic_charge=True,charge_square_pairs=True,full_charge_indicators=True,full_signed_charge=True)
    # With one sparse shape the old bound is58 and the enlarged bound is59.
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**flags)
    assert _replay(c,**flags,hopping_telescope=True)['mixture_sources']==59
    flags['full_signed_charge']=False
    with pytest.raises(ValueError,match='requires full signed'):_replay(c,**flags,hopping_telescope=True)
    flags['full_signed_charge']=True
    with pytest.raises(ValueError,match='requires full signed'):_replay(c,**flags,hopping_telescope=1)
