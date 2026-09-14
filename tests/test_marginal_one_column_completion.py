"""Exact feasibility and refusal tests for physical-column completion."""
from pathlib import Path
from fractions import Fraction as F
import sys
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/discovery'))
from one_column_physical_completion import complete


def test_opposite_physical_moment_completes_normalized_mixture():
    index, weights, _ = complete([[1],[1]],[1,0],[(7,[1,-1])])
    assert index == 7 and weights == [F(1,2),F(1,2)]


def test_two_moments_completed_simultaneously_with_positive_mass():
    _, weights, _ = complete([[1,1],[1,0],[0,1]],[1,0,0],[(3,[1,-1,-1])])
    assert weights == [F(1,3)]*3


def test_exact_completion_resolves_residual_below_float_precision():
    tiny = F(1,10**30)
    _, weights, _ = complete([[1],[1]],[1,1+tiny],[(7,[1,2])])
    assert weights == [1-tiny,tiny]


def test_completion_refuses_negative_existing_or_new_mass():
    with pytest.raises(ValueError,match='No nonnegative physical completion'):
        complete([[1],[1]],[1,2],[(0,[1,F(3,2)])])
    with pytest.raises(ValueError,match='No nonnegative physical completion'):
        complete([[1],[1]],[1,0],[(0,[1,2])])


def test_dependent_original_columns_and_wrong_codimension_are_refused():
    with pytest.raises(ValueError,match='Dependent'):
        complete([[1,2],[1,2],[1,2]],[1,0,0],[(0,[1,-1,0])])
    with pytest.raises(ValueError,match='One missing column'):
        complete([[1],[1],[1]],[1,0,0],[(0,[1,-1,0])])
