"""Independent contractions and fermionic checks for compact overlap telescopes."""
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/discovery'))
from full_overlap_telescope import build, embed, permutation, shift6, moment, hole5
from full_overlap_density import reconstruct, quadratic
from experiments.marginal_transfer_verify import apply_word


def test_fermionic_cyclic_translation_matches_creation_operator_construction():
    for state in range(4096):
        target, sign = 0, 1
        # Canonical Fock states are products of ascending creation modes, acting right to left.
        for mode in reversed(range(12)):
            if state >> mode & 1:
                target, phase = apply_word(((1, (mode+2) % 12),), target)
                sign *= phase
        assert shift6(state) == (target, sign)


def test_contiguous_embeddings_match_kronecker_products():
    small = {(3, 5): 2, (5, 3): 2, (6, 6): -3}
    dense = np.zeros((1024, 1024), dtype=np.int8)
    for (row, column), value in small.items():
        dense[row, column] = value
    for left in (True, False):
        expected = np.kron(np.eye(4, dtype=np.int8), dense) if left else np.kron(dense, np.eye(4, dtype=np.int8))
        rows, columns = np.nonzero(expected)
        assert embed(small, left) == {(int(r), int(c)): int(expected[r, c]) for r, c in zip(rows, columns)}


@pytest.mark.parametrize('vector', [{358: 1, 601: -1}, {346: 1, 613: 1}])
def test_direct_moment_matches_independent_symmetry_averaged_partial_traces(vector):
    denominator, odd, telescope = build(vector)
    mixture = [{'weight': '2/7', 'vector': {'1366': 2, '1621': -3}},
               {'weight': '5/7', 'vector': {'346': 1, '613': 4}}]
    scale, left, right, difference, _ = reconstruct(mixture)
    expected = F(quadratic(difference, vector), scale*sum(a*a for a in vector.values()))
    assert moment(mixture, telescope, denominator) == expected
    assert telescope and odd
    # Numerical eigenvalues are only an independent test; the receipt uses the projector proof.
    states = sorted({s for key in telescope for s in key}); indices = {s: i for i, s in enumerate(states)}
    matrix = np.zeros((len(states), len(states)))
    for (r, c), a in telescope.items():
        matrix[indices[r], indices[c]] = a/denominator
    assert np.linalg.norm(matrix, 2) <= 1+1e-13


def test_odd_site_particle_hole_is_projective_not_vector_involution():
    for state in range(1024):
        target, sign = hole5(state); restored, phase = hole5(target)
        assert restored == state and sign*phase == -1


@pytest.mark.parametrize('vector', [{}, {1024: 1}, {0: 1, 1: 1}, {0: 0},
                                   {0: 1, 1: 2, 2: 3}, {0: True}, {0: 10**9+1}])
def test_invalid_sources_refused(vector):
    with pytest.raises(ValueError):
        build(vector)
