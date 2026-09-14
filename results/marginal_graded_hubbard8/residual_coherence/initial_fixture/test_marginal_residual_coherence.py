"""Independent partial traces and refusal checks for pure-coherence obstructions."""
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'results/marginal_graded_hubbard8/discovery'))
from residual_coherence_obstruction import pure_coherence, functional, verify_annihilators, discover_functionals, eligible_entry, OLD_MODULES
import importlib

FUNCTIONALS = (((413, 1433, F(-16)),), ((1382, 1433, F(8)),))
from full_overlap_telescope import build, moment, embed
from full_overlap_density import reconstruct

VECTORS = ({358: 1, 409: 1}, {103: 1, 358: 1})


@pytest.mark.parametrize('vector', VECTORS)
def test_pure_moment_matches_independent_offdiagonal_partial_trace(vector):
    mixture = [{'weight': '2/7', 'vector': {str(s): 2 if i == 0 else -3 for i, s in enumerate(vector)}},
               {'weight': '5/7', 'vector': {'346': 1, '409': 4}}]
    denominator, y, telescope, bound = pure_coherence(vector)
    scale, left, right, difference, images = reconstruct(mixture)
    (r, a), (c, b) = sorted(vector.items())
    expected = F(2 * a * b * difference.get((r, c), 0), scale * sum(v*v for v in vector.values()))
    assert moment(mixture, telescope, denominator) == expected
    # Removing the diagonal is a matrix operation, not permission to ignore its expectation.
    full_denominator, _, full = build(vector)
    diagonal = {key: value for key, value in full.items() if key[0] == key[1]}
    assert moment(mixture, full, full_denominator) == expected + moment(mixture, diagonal, full_denominator)
    assert not any(r == c for r, c in telescope)
    assert bound == (F(1, 4) if 409 in vector else F(1, 8))


@pytest.mark.parametrize('vector', VECTORS)
def test_numerical_norm_independently_checks_exact_row_bound(vector):
    denominator, y, telescope, bound = pure_coherence(vector)
    for matrix, limit in ((y, bound), (telescope, 2 * bound)):
        states = sorted({s for key in matrix for s in key})
        indices = {s: i for i, s in enumerate(states)}
        dense = np.zeros((len(states), len(states)))
        for (r, c), value in matrix.items():
            dense[indices[r], indices[c]] = value / denominator
        assert np.max(abs(np.linalg.eigvalsh(dense))) <= float(limit) + 1e-13


@pytest.mark.parametrize('vector', [{0: 1}, {346: 1}, {0: 1, 1: 1}, {346: True, 409: 1}])
def test_empty_coherence_and_invalid_physical_source_refused(vector):
    with pytest.raises(ValueError):
        pure_coherence(vector)


def test_exact_minor_and_annihilator_refusals():
    matrices = []
    for vector in VECTORS:
        denominator, _, t, _ = pure_coherence(vector)
        matrices.append({key: F(value, denominator) for key, value in t.items()})
    minor, determinant = verify_annihilators([], matrices, FUNCTIONALS)
    assert minor == [[F(0), F(1)], [F(1), F(0)]]
    assert determinant == -1
    with pytest.raises(ValueError, match='independent'):
        verify_annihilators([], [matrices[0], matrices[0]], FUNCTIONALS)
    fake = [{} for _ in range(4096)]
    r, c, sign = FUNCTIONALS[0][0]
    fake[c][r] = F(1)
    with pytest.raises(ValueError, match='annihilate'):
        verify_annihilators([('injected-old-support', fake)], matrices, FUNCTIONALS)
    with pytest.raises(ValueError, match='4096'):
        verify_annihilators([('incomplete-action', fake[:-1])], matrices, FUNCTIONALS)


def test_discovery_against_all_existing_actions_and_dependent_refusal():
    old = []
    for name in OLD_MODULES:
        module = importlib.import_module('experiments.marginal_' + name)
        for label in getattr(module, 'LABELS', ('single',)):
            old.append((name + ':' + label, module.actions({label: 1}) if hasattr(module, 'LABELS') else module.actions()))
    matrices = []
    for vector in VECTORS:
        denominator, _, t, _ = pure_coherence(vector)
        matrices.append({key:F(value,denominator) for key,value in t.items()})
    found = discover_functionals(old, matrices)
    assert found == FUNCTIONALS
    assert verify_annihilators(old, matrices, found)[1] == -1
    with pytest.raises(ValueError, match='independent'):
        discover_functionals(old, [matrices[0],matrices[0]])
    with pytest.raises(ValueError, match='Bounded'):
        discover_functionals(old[:-1], matrices)


def test_two_bit_long_hop_excludes_actual_nearest_neighbor_physical_actions():
    from experiments.marginal_local_hubbard_block import _actions
    action = _actions(6, F(10,3), 1, [F(i+1,7) for i in range(6)],
                      [F(2*i+1,11) for i in range(5)], F(1,2), [F(i+1,13) for i in range(5)])
    for terms in FUNCTIONALS:
        assert all(eligible_entry(r,c) for r,c,a in terms)
        assert sum(a*action[c].get(r,0) for r,c,a in terms) == 0
    # Same spin sector is necessary, but nearest-neighbor support remains ineligible.
    assert not eligible_entry(1370, 1370)
    assert not eligible_entry(0, 1)
