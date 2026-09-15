"""Exact coherent-projector integration, physical coverage and refusal gates."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import subprocess
import pytest
from experiments.marginal_coherent_projector_telescope import LABELS,coefficients,actions
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_range_two_family_limit import replay as family_replay
from experiments.marginal_diagonal_family_limit import _replay
from experiments.marginal_hopping_telescope import projected_matrix
from experiments.marginal_local_hubbard_block import sector_matrices
from test_marginal_three_spectator_hopping import three_toy
from test_marginal_spectator_hopping import family
from test_marginal_pair_family import flags

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'results/marginal_graded_hubbard8/discovery'))
from full_overlap_telescope import build as independent


def toy():
    c=three_toy();c['kind']='hubbard_projector_extension_v17'
    c['coherent_projector']={label:'1/1000' for label in LABELS}
    return c


def test_production_operators_match_independent_exact_overlap_construction():
    for label,vector in LABELS.items():
        denominator,_,matrix=independent(vector)
        actual=actions({label:1})
        assert {(r,c):a for c,image in enumerate(actual) for r,a in image.items()}=={key:F(a,denominator) for key,a in matrix.items()}


def test_mixed_projection_covers_every_physical_basis_image():
    data=actions(dict(zip(LABELS,(F(2,7),F(-3,11)))))
    for _,_,cols in sector_matrices(6,0,0).values():
        matrix=projected_matrix(cols,data)
        for j,col in enumerate(cols):
            actual={}
            for s,a in col.items():
                for t,b in data[s].items():actual[t]=actual.get(t,0)+a*b
            rebuilt={s:F(matrix[i][j]*a,sum(v*v for v in basis.values()))
                     for i,basis in enumerate(cols) for s,a in basis.items()}
            assert {s:a for s,a in actual.items() if a}=={s:a for s,a in rebuilt.items() if a}


def test_energy_accepts_safe_case_and_refuses_negative_psd():
    c=toy();r=replay(c)
    assert r['accepted'] and r['coherent_projector']==c['coherent_projector']
    assert r['local_sum_dimensions']==4096 and len(r['local_sectors'])==94
    assert r['open_lower_density']==r['periodic_lower_density']=='-9/5'
    c['coherent_projector']={next(iter(LABELS)):'100'}
    with pytest.raises(ValueError):replay(c)


@pytest.mark.parametrize('version',[1,2,4,5,6,7,8,9,10,11,12,13,14,15,16])
def test_old_versions_cannot_ignore_new_operator(version):
    c=toy();c['kind']=f'hubbard_projector_extension_v{version}'
    with pytest.raises(ValueError,match='requires v17'):replay(c)


@pytest.mark.parametrize('value',[None,{}, {'unknown':1},{next(iter(LABELS)):0},
                                {next(iter(LABELS)):True},{next(iter(LABELS)):0.1},
                                {next(iter(LABELS)):'1/1000001'}])
def test_invalid_coefficients_refused(value):
    with pytest.raises(ValueError):coefficients(value)


def test_stationary_family_and_fixed_coefficient_refusal():
    c=family();c['kind']='joint_diagonal_range2_family_limit_v13';r=family_replay(c)
    assert r['coherent_projector_moments']==dict.fromkeys(LABELS,'0')
    assert len(r['three_spectator_hopping_moments'])==18 and r['periodic_family_upper']=='13/5'
    c['coherent_projector']={next(iter(LABELS)):'1/1000'}
    with pytest.raises(ValueError,match='constrain moments'):family_replay(c)


@pytest.mark.parametrize('case',['W_zero','W_plus_1'])
def test_old_nonstationary_mixtures_cannot_be_relabelled(case):
    c=json.loads((ROOT/f'results/marginal_graded_hubbard8/three_spectator/{case}/final/range_two_family_limit_certificate.json').read_text())
    c['kind']='joint_diagonal_range2_family_limit_v13'
    with pytest.raises(ValueError,match='coherent-projector expectation'):family_replay(c)


def test_bounded_source_count_and_hierarchy():
    c=family();c.pop('W');c.pop('range_two_density_profile');c['kind']='joint_diagonal_family_limit_v1'
    f=flags();f.update(pair_transfer=True,two_spectator_hopping=True,three_spectator_hopping=True)
    c['mixture']=[{'weight':'1/131','vector':{0:1}}]*131
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**f)
    assert _replay(c,**f,coherent_projector=True)['mixture_sources']==131
    c['mixture']=[{'weight':'1/132','vector':{0:1}}]*132
    with pytest.raises(ValueError,match='constraint count'):_replay(c,**f,coherent_projector=True)
    f['three_spectator_hopping']=False
    with pytest.raises(ValueError,match='requires three-spectator'):_replay(c,**f,coherent_projector=True)
    f['three_spectator_hopping']=True
    with pytest.raises(ValueError,match='requires three-spectator'):_replay(c,**f,coherent_projector=1)


def test_matching_driver_refuses_old_family_for_new_energy(tmp_path):
    p=ROOT/'results/marginal_graded_hubbard8/three_spectator/W_zero/final'
    (tmp_path/'diagonal_family_limit_proposal.json').write_bytes((p/'diagonal_family_limit_proposal.json').read_bytes())
    c=json.loads((p/'profile_joint_r1_2_certificate.json').read_text());c['kind']='hubbard_projector_extension_v17'
    c['coherent_projector']={next(iter(LABELS)):'1/1000'}
    (tmp_path/'profile_joint_r1_2_certificate.json').write_text(json.dumps(c))
    result=subprocess.run([sys.executable,str(ROOT/'results/marginal_graded_hubbard8/discovery/range_two_family_limit_replay.py'),str(tmp_path)],capture_output=True,text=True,timeout=90)
    assert result.returncode!=0 and 'Coherent-projector energy requires a matching enlarged family cap' in result.stderr
    assert not (tmp_path/'range_two_family_limit_replay.json').exists()
