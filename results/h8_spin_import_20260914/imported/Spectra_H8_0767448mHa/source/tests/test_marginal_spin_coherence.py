"""Independent partial traces and refusal checks for pure-coherence obstructions."""
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'results/marginal_graded_hubbard8/discovery'))
from spin_coherence_obstruction import pure_coherence, functional, verify_annihilators, FUNCTIONALS
from full_overlap_telescope import build, moment, embed
from full_overlap_density import reconstruct

VECTORS = ({346: 1, 409: 1}, {314: 1, 614: 1})


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
    assert bound == F(1, 4)


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
    minor, determinant = verify_annihilators([], matrices)
    assert minor == [[F(1, 8), F(0)], [F(0), F(1, 16)]]
    assert determinant == F(1, 128)
    with pytest.raises(ValueError, match='independent'):
        verify_annihilators([], [matrices[0], matrices[0]])
    fake = [{} for _ in range(4096)]
    r, c, sign = FUNCTIONALS[0][0]
    fake[c][r] = F(1)
    with pytest.raises(ValueError, match='annihilate'):
        verify_annihilators([('injected-old-support', fake)], matrices)
    with pytest.raises(ValueError, match='4096'):
        verify_annihilators([('incomplete-action', fake[:-1])], matrices)
