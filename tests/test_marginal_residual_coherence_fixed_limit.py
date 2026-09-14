"""Fixed-recipe cap proof, actual amplitudes and refusal paths."""
from pathlib import Path
from fractions import Fraction as F
import json
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'results/marginal_graded_hubbard8/discovery'))
from residual_coherence_fixed_limit_replay import replay, fixed_expectation


def seed(case='W_zero'):
    selection = 'refined' if case == 'W_zero' else 'scaled'
    return json.loads((ROOT / f'results/marginal_graded_hubbard8/pure_coherence/{case}/{selection}/profile_joint_r1_2_certificate.json').read_text())


def proposal(mixture=None):
    return {'kind': 'residual_coherence_fixed_family_proposal_v1',
            'mixture': [{'weight': '1', 'vector': {'0': 1}}] if mixture is None else mixture}


@pytest.mark.parametrize('case', ['W_zero', 'W_plus_1'])
def test_vacuum_expectation_matches_direct_occupation_energy(case):
    c = seed(case)
    result = replay(c, proposal())
    expected = sum(map(F, c['local_window']['onsite_profile']))/2
    expected += sum(map(F, c['local_window']['density_profile']))
    expected += sum(map(F, c['local_window']['range_two_density_profile']))
    assert result['exact_source_local_expectations'] == [str(expected)]
    assert result['new_coherence_moments'] == ['0', '0']
    assert result['new_coefficients_unrestricted_real'] and result['all_old_fields_frozen']
    scaled = replay(c, proposal([{'weight': '1', 'vector': {'0': -7}}]))
    assert scaled == result


@pytest.mark.parametrize('mixture', [[], [{'weight': '1/2', 'vector': {'0': 1}}],
                                   [{'weight': '-1', 'vector': {'0': 1}}],
                                   [{'weight': '1', 'vector': {'0': 1, '1': 1}}],
                                   [{'weight': '1', 'vector': {'358': 1, '409': 1}}],
                                   [{'weight': '1/4', 'vector': {'0': 1}}]*4])
def test_invalid_or_nonstationary_mixture_refused(mixture):
    with pytest.raises(ValueError):
        replay(seed(), proposal(mixture))


def test_wrong_family_scope_refused():
    p = proposal()
    p['kind'] = 'joint_diagonal_range2_family_limit_v14'
    with pytest.raises(ValueError, match='frozen-recipe'):
        replay(seed(), p)
    c = seed()
    c['target']['U'] = '5'
    with pytest.raises(ValueError, match='Matched'):
        fixed_expectation(c)


def test_frozen_expectation_includes_preceding_pure_coherence_term():
    from copy import deepcopy
    c = seed()
    without = deepcopy(c)
    without['pure_coherence'] = {label:'0' for label in c['pure_coherence']}
    vector = {1370:1,1433:1}
    difference = fixed_expectation(c)(vector,2)-fixed_expectation(without)(vector,2)
    assert difference == F(c['pure_coherence']['346,409,1'])/8
    assert difference


def test_preceding_recipe_version_cannot_be_relabelled():
    c = seed()
    c['kind'] = 'hubbard_projector_extension_v18'
    with pytest.raises(ValueError, match='ENERGYv19'):
        fixed_expectation(c)
