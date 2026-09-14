"""Production integration of the two exact residual-coherence directions."""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import subprocess
import pytest
from experiments.marginal_residual_coherence import LABELS, build, coefficients, actions
from experiments.marginal_projector_extendibility import replay
from experiments.marginal_range_two_family_limit import replay as family_replay
from experiments.marginal_diagonal_family_limit import _replay
from test_marginal_pure_coherence import toy as word_toy
from test_marginal_spectator_hopping import family
from test_marginal_pair_family import flags

ROOT = Path(__file__).resolve().parents[1]


def toy():
    c = word_toy()
    c['kind'] = 'hubbard_projector_extension_v20'
    c['residual_coherence'] = {key: '1/1000' for key in LABELS}
    return c


@pytest.mark.parametrize('case,label', list(zip(('W_zero', 'W_plus_1'), LABELS)))
def test_matches_independently_accepted_sparse_operator(case, label):
    source = json.loads((ROOT / f'results/marginal_graded_hubbard8/residual_coherence/{case}/telescope_replay.json').read_text())
    denominator, y, t = build(label)
    assert denominator == source['matrix_denominator']
    assert [[r, c, a] for (r, c), a in sorted(y.items())] == source['five_site_numerator']
    assert [[r, c, a] for (r, c), a in sorted(t.items())] == source['six_site_numerator']
    action = actions({label: F(1)})
    assert {(r, c): a for c, image in enumerate(action) for r, a in image.items()} == {key: F(value, denominator) for key, value in t.items()}
    assert all(r != c for r, c in t)


def test_actual_energy_acceptance_and_negative_psd_refusal():
    c = toy()
    r = replay(c)
    assert r['accepted'] and r['residual_coherence'] == c['residual_coherence']
    assert r['local_sum_dimensions'] == 4096 and len(r['local_sectors']) == 94
    assert r['periodic_lower_density'] == r['open_lower_density'] == '-9/5'
    c['residual_coherence'] = {next(iter(LABELS)): '1000'}
    with pytest.raises(ValueError):
        replay(c)


@pytest.mark.parametrize('version', range(1, 20))
def test_all_older_versions_refuse_new_field(version):
    c = toy()
    c['kind'] = f'hubbard_projector_extension_v{version}'
    with pytest.raises(ValueError, match='requires v20'):
        replay(c)


@pytest.mark.parametrize('value', [None, {}, {'bad': 1}, {next(iter(LABELS)): 0},
                                {next(iter(LABELS)): True}, {next(iter(LABELS)): .1},
                                {next(iter(LABELS)): '1/1000001'}])
def test_exact_coefficient_refusals(value):
    with pytest.raises(ValueError):
        coefficients(value)


def test_exact_family_moments_and_fixed_field_refusal():
    c = family()
    c['kind'] = 'joint_diagonal_range2_family_limit_v16'
    r = family_replay(c)
    assert r['residual_coherence_moments'] == dict.fromkeys(LABELS, '0')
    assert len(r['spin_word_moments']) == 120 and r['periodic_family_upper'] == '13/5'
    c['residual_coherence'] = {next(iter(LABELS)): '1/1000'}
    with pytest.raises(ValueError, match='constrain moments'):
        family_replay(c)


@pytest.mark.parametrize('case,selection', [('W_zero', 'refined'), ('W_plus_1', 'scaled')])
def test_previous_nonstationary_mixture_cannot_be_relabelled(case, selection):
    c = json.loads((ROOT / f'results/marginal_graded_hubbard8/pure_coherence/{case}/{selection}/range_two_family_limit_certificate.json').read_text())
    c['kind'] = 'joint_diagonal_range2_family_limit_v16'
    with pytest.raises(ValueError, match='residual-coherence expectation'):
        family_replay(c)


def test_preserved_old_source_cap_and_new_hierarchy():
    c = family()
    c.pop('W'); c.pop('range_two_density_profile')
    c['kind'] = 'joint_diagonal_family_limit_v1'
    f = flags()
    f.update(pair_transfer=True, two_spectator_hopping=True, three_spectator_hopping=True,
             coherent_projector=True, full_spin_word=True, pure_coherence=True)
    c['mixture'] = [{'weight': '1/203', 'vector': {0: 1}}]*203
    with pytest.raises(ValueError, match='constraint count'):
        _replay(c, **f)
    assert _replay(c, **f, residual_coherence=True)['mixture_sources'] == 203
    c['mixture'] = [{'weight': '1/204', 'vector': {0: 1}}]*204
    with pytest.raises(ValueError, match='constraint count'):
        _replay(c, **f, residual_coherence=True)
    f['pure_coherence'] = False
    with pytest.raises(ValueError, match='requires pure-coherence hierarchy'):
        _replay(c, **f, residual_coherence=True)
    with pytest.raises(ValueError, match='requires pure-coherence hierarchy'):
        _replay(c, **f, residual_coherence=1)


def test_old_family_driver_refuses_new_energy(tmp_path):
    source = ROOT / 'results/marginal_graded_hubbard8/pure_coherence/W_zero/refined'
    (tmp_path / 'diagonal_family_limit_proposal.json').write_bytes((source / 'diagonal_family_limit_proposal.json').read_bytes())
    c = json.loads((source / 'profile_joint_r1_2_certificate.json').read_text())
    c['kind'] = 'hubbard_projector_extension_v20'
    c['residual_coherence'] = {next(iter(LABELS)): '1/1000'}
    (tmp_path / 'profile_joint_r1_2_certificate.json').write_text(json.dumps(c))
    result = subprocess.run([sys.executable, str(ROOT / 'results/marginal_graded_hubbard8/discovery/range_two_family_limit_replay.py'), str(tmp_path)], capture_output=True, text=True, timeout=180)
    assert result.returncode != 0 and 'Residual-coherence energy requires a matching enlarged family cap' in result.stderr
    assert not (tmp_path / 'range_two_family_limit_replay.json').exists()
