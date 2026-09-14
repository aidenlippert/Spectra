"""Fixed-recipe cap proof, actual amplitudes and refusal paths."""
from pathlib import Path
from fractions import Fraction as F
import json
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'results/marginal_graded_hubbard8/discovery'))
from spin_coherence_fixed_limit_replay import replay, fixed_expectation


def seed(case='W_zero'):
    selection = 'scaled' if case == 'W_zero' else 'final'
    return json.loads((ROOT / f'results/marginal_graded_hubbard8/spin_word/{case}/{selection}/profile_joint_r1_2_certificate.json').read_text())


def proposal(mixture=None):
    return {'kind': 'spin_coherence_fixed_family_proposal_v1',
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
                                   [{'weight': '1', 'vector': {'346': 1, '409': 1}}],
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
