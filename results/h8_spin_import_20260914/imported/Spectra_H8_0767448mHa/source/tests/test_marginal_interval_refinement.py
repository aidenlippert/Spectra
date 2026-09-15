"""Refinement comparisons must preserve the family and both exact endpoints."""
import importlib.util
from pathlib import Path
import pytest

path = Path(__file__).resolve().parents[1] / 'results/marginal_graded_hubbard8/discovery/spectator_interval_refinement.py'
spec = importlib.util.spec_from_file_location('interval_refinement', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def receipt(lower='1', upper='3', gap='2'):
    return dict(accepted=True, spectator_hopping=True,
                accepted_periodic_lower=lower, periodic_family_upper=upper, family_gap=gap)


def certificate():
    return dict.fromkeys(module.FIXED, 'fixed')


def test_exact_nested_refinement():
    result = module.compare(receipt(), receipt('3/2', '5/2', '1'), certificate(), certificate())
    assert result['gap_reduction_factor'] == '2'


@pytest.mark.parametrize('new', [receipt('0', '1', '1'), receipt('3', '4', '1'),
    receipt(), receipt('2', '1', '-1'), receipt('3/2', '5/2', '1/2'),
    dict(receipt('2', '3', '1'), accepted=False)])
def test_refuse_invalid_or_non_nested(new):
    with pytest.raises(ValueError):
        module.compare(receipt(), new, certificate(), certificate())


@pytest.mark.parametrize('key', module.FIXED)
def test_refuse_changed_family(key):
    new_certificate = certificate()
    new_certificate[key] = 'changed'
    with pytest.raises(ValueError):
        module.compare(receipt(), receipt('2', '3', '1'), certificate(), new_certificate)
