import json
from pathlib import Path
import pytest

from experiments.marginal_diagonal_family_limit import replay


def certificate():
    return {'kind': 'joint_diagonal_family_limit_v1', 'half_vector': {0x999: 1, 0x666: 1},
            'charged_vector': {62: 1, 3008: 1}, 'ratio': '1/2', 'theta_half': '1/2',
            'theta_joint': '3/4', 'diagonal_shapes': [{102: 1, 612: -1}],
            'mixture': [{'weight': '1/3', 'vector': {0: 7}}, {'weight': '2/3', 'vector': {4095: 11}}]}


def test_vacuum_full_mixture_caps_the_whole_diagonal_span():
    result = replay(certificate())
    assert result['periodic_family_upper'] == '5/2'
    assert result['diagonal_moments'] == ['0']
    assert result['profile_gradients'] == ['0']*6


def test_old_accepted_dual_is_refused_for_a_violated_diagonal():
    root = Path(__file__).resolve().parents[1]
    c = json.loads((root/'results/marginal_graded_hubbard8/joint_projector/signed_density/family_limit_certificate.json').read_text())
    c.update(kind='joint_diagonal_family_limit_v1', diagonal_shapes=[{358: 1, 613: -1}])
    with pytest.raises(ValueError, match='diagonal expectation'):
        replay(c)


@pytest.mark.parametrize('field,value', [('ratio', '0'), ('theta_half', '2'),
                                        ('diagonal_shapes', []), ('diagonal_shapes', [{102: 1}]),
                                        ('diagonal_shapes', [{102: 1, 612: -1}]*10),
                                        ('proposed_periodic_family_upper', '0'),
                                        ('kind', 'joint_diagonal_family_limit_proposal_v1')])
def test_malformed_claims_refused(field, value):
    c = certificate()
    c[field] = value
    with pytest.raises(ValueError):
        replay(c)


def test_mixture_trace_positivity_size_and_fidelity_refusals():
    for mixture in ([{'weight': '1/2', 'vector': {0: 1}}],
                    [{'weight': '-1', 'vector': {0: 1}}, {'weight': '2', 'vector': {4095: 1}}],
                    [{'weight': '1', 'vector': {0x999: 1, 0x666: 1}}],
                    [{'weight': '1/11', 'vector': {0: 1}}]*11):
        c = certificate()
        c['mixture'] = mixture
        with pytest.raises(ValueError):
            replay(c)


def test_shape_union_cannot_bypass_energy_entry_cap():
    from experiments.marginal_local_hubbard_block import _reflection
    pairs = [(s, _reflection(s, 5)[0]) for s in range(1024) if s < _reflection(s, 5)[0]][:33]
    shapes = []
    for start in range(0, 33, 11):
        shape = {}
        for s, reflected in pairs[start:start+11]:
            shape.update({s: 1, reflected: -1})
        shapes.append(shape)
    c = certificate()
    c['diagonal_shapes'] = shapes
    with pytest.raises(ValueError, match='64-entry'):
        replay(c)
