"""Independent exact checks for the bounded discovery solver."""
from fractions import Fraction as F
from pathlib import Path
import sys
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/marginal_graded_hubbard8/discovery'))
from spin_word_fraction_free import solve
from two_spectator_family_resume import solve_basis


def test_fraction_free_matches_independent_fraction_elimination():
    matrix=[[F(1,i+j+1) for j in range(9)] for i in range(11)]
    rhs=[sum(v*F(1,j+1) for j,v in enumerate(row)) for row in matrix]
    assert solve(matrix,rhs)==solve_basis(matrix,rhs)==[F(1,j+1) for j in range(9)]


def test_full_row_bound_and_rational_column_scalings():
    matrix=[[F(1,3),F(2,5)],[1,-1]]+[[F(1,3),F(2,5)]]*196
    assert solve(matrix,[1,0]+[1]*196)==[F(15,11),F(15,11)]


def test_inconsistent_dependent_and_over_budget_systems_are_refused():
    with pytest.raises(ValueError,match='inconsistent'):solve([[1],[1]],[1,2])
    with pytest.raises(ValueError,match='Dependent'):solve([[1,2],[2,4]],[1,2])
    with pytest.raises(ValueError,match='198-row'):solve([[1]]*199,[1]*199)
    with pytest.raises(ValueError,match='4096-bit'):solve([[2**4097]],[1])
