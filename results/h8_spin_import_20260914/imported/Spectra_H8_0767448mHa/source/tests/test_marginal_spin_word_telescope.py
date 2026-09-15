"""Complete spin-word basis, production PSD integration and preserved gates."""
from fractions import Fraction as F
from pathlib import Path
import json,subprocess,sys
import pytest
from experiments.marginal_spin_word_telescope import PATTERNS,ORBIT,coefficients,five_site_value,local_value,spin_flip
from experiments.marginal_local_hubbard_block import _reflection
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_range_two_family_limit import replay as family_replay
from experiments.marginal_diagonal_family_limit import _replay
from experiments.marginal_signed_charge_telescope import PATTERNS as CHARGE,five_site_value as charge_value
from test_marginal_coherent_projector import toy as coherent_toy
from test_marginal_spectator_hopping import family
from test_marginal_pair_family import flags

ROOT=Path(__file__).resolve().parents[1]

def toy():
    c=coherent_toy();c['kind']='hubbard_projector_extension_v18'
    c['spin_word_telescope']={key:'1/1000' for key in PATTERNS}
    return c

def test_basis_reconstructs_independent_full_group_projection():
    f=[(s*s+13*s+7)%31-15 for s in range(1024)]
    projected={}
    for s in range(1024):
        r=_reflection(s,5)[0]
        projected[s]=F(sum(f[t] for t in (s,1023^s,spin_flip(s),1023^spin_flip(s)))-sum(f[t] for t in (r,1023^r,spin_flip(r),1023^spin_flip(r))),8)
    terms={key:projected[int(key)] for key in PATTERNS}
    assert len(PATTERNS)==120 and sum(len(v) for v in PATTERNS.values())==len(ORBIT)
    assert all(five_site_value(s,terms)==projected[s] for s in range(1024))

def test_every_old_charge_direction_and_selected_sparse_shape_is_in_full_space():
    vectors=[{s:charge_value(s,{key:F(1)}) for s in range(1024)} for key in CHARGE]
    c=json.loads((ROOT/'results/marginal_graded_hubbard8/coherent_projector/W_zero/final/range_two_family_limit_certificate.json').read_text())
    vectors += [{int(s):F(v) for s,v in shape.items()} for shape in c['diagonal_shapes']]
    for vector in vectors:
        terms={key:vector.get(int(key),F(0)) for key in PATTERNS}
        assert all(five_site_value(s,terms)==vector.get(s,F(0)) for s in range(1024))

def test_every_fock_state_has_exact_periodic_telescoping_and_reflection():
    terms={key:F(i%13-6,17) for i,key in enumerate(PATTERNS)}
    for s in range(4096):
        current=s;total=F(0)
        for _ in range(6):
            total+=local_value(current,terms);current=((current<<2)&4095)|(current>>10)
        assert current==s and total==0
        assert local_value(_reflection(s,6)[0],terms)==local_value(s,terms)

def test_safe_energy_and_negative_psd_refusal():
    c=toy();r=replay(c)
    assert r['accepted'] and r['spin_word_telescope']==c['spin_word_telescope']
    assert r['local_sum_dimensions']==4096 and len(r['local_sectors'])==94
    assert r['open_lower_density']==r['periodic_lower_density']=='-9/5'
    c['spin_word_telescope']={next(iter(PATTERNS)):'100'}
    with pytest.raises(ValueError):replay(c)

@pytest.mark.parametrize('version',[1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
def test_old_energy_versions_cannot_ignore_spin_words(version):
    c=toy();c['kind']=f'hubbard_projector_extension_v{version}'
    with pytest.raises(ValueError,match='requires v18'):replay(c)

@pytest.mark.parametrize('value',[None,{}, {'unknown':1},{next(iter(PATTERNS)):0},
                                {next(iter(PATTERNS)):True},{next(iter(PATTERNS)):0.1},
                                {next(iter(PATTERNS)):'1/1000001'}])
def test_invalid_exact_coefficients_refused(value):
    with pytest.raises(ValueError):coefficients(value)

def test_stationary_family_all120_moments_and_fixed_coefficient_refusal():
    c=family();c['kind']='joint_diagonal_range2_family_limit_v14';r=family_replay(c)
    assert r['spin_word_moments']==dict.fromkeys(PATTERNS,'0')
    assert len(r['coherent_projector_moments'])==2 and r['periodic_family_upper']=='13/5'
    c['spin_word_telescope']={next(iter(PATTERNS)):'1/1000'}
    with pytest.raises(ValueError,match='constrain moments'):family_replay(c)

@pytest.mark.parametrize('case',['W_zero','W_plus_1'])
def test_violated_previous_mixture_cannot_be_relabelled(case):
    c=json.loads((ROOT/f'results/marginal_graded_hubbard8/coherent_projector/{case}/final/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v14'
    with pytest.raises(ValueError,match='spin-word expectation'):family_replay(c)

def test_new_source_cap_and_required_hierarchy():
    c=family();c.pop('W');c.pop('range_two_density_profile');c['kind']='joint_diagonal_family_limit_v1'
    f=flags();f.update(pair_transfer=True,two_spectator_hopping=True,three_spectator_hopping=True,coherent_projector=True)
    c['mixture']=[{'weight':'1/199','vector':{0:1}}]*199
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**f)
    assert _replay(c,**f,full_spin_word=True)['mixture_sources']==199
    c['mixture']=[{'weight':'1/200','vector':{0:1}}]*200
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**f,full_spin_word=True)
    f['coherent_projector']=False
    with pytest.raises(ValueError,match='requires coherent-projector'):_replay(c,**f,full_spin_word=True)
    f['coherent_projector']=True
    with pytest.raises(ValueError,match='requires coherent-projector'):_replay(c,**f,full_spin_word=1)

def test_old_family_driver_cannot_accept_new_energy(tmp_path):
    p=ROOT/'results/marginal_graded_hubbard8/coherent_projector/W_zero/final'
    (tmp_path/'diagonal_family_limit_proposal.json').write_bytes((p/'diagonal_family_limit_proposal.json').read_bytes())
    c=json.loads((p/'profile_joint_r1_2_certificate.json').read_text());c['kind']='hubbard_projector_extension_v18'
    c['spin_word_telescope']={next(iter(PATTERNS)):'1/1000'}
    (tmp_path/'profile_joint_r1_2_certificate.json').write_text(json.dumps(c))
    r=subprocess.run([sys.executable,str(ROOT/'results/marginal_graded_hubbard8/discovery/range_two_family_limit_replay.py'),str(tmp_path)],capture_output=True,text=True,timeout=120)
    assert r.returncode!=0 and 'Full spin-word energy requires a matching enlarged family cap' in r.stderr
    assert not (tmp_path/'range_two_family_limit_replay.json').exists()
