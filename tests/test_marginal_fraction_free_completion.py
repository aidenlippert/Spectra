import importlib.util
from fractions import Fraction as F
from pathlib import Path

import pytest


path = Path(__file__).resolve().parents[1] / 'results/marginal_graded_hubbard8/discovery/fraction_free_completion.py'
spec = importlib.util.spec_from_file_location('fraction_free_completion', path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
reduce_basis = module.reduce_basis


@pytest.mark.parametrize('columns,rhs', [
    ([[2, 1], [1, 3]], [5, 7]),
    ([[1, 2], [3, 1]], [5, 7]),
    ([[0, 2], [3, 1]], [4, 9]),
    ([[1, 2, 3]], [4, 8, 12]),
    ([[0, 2, 3], [1, 0, 4]], [5, 8, 32]),
    ([[0, 2, 3], [1, 0, 4]], [5, 8, 33]),
])
def test_reduction_reconstructs_original_equations(columns, rhs):
    # Solve against every coordinate rhs, then verify the induced row map on
    # original columns and rhs. This also checks rectangular residual rows.
    m, n = len(rhs), len(columns)
    identity = [[int(i == j) for i in range(m)] for j in range(m)]
    d, b, transform = reduce_basis(columns, rhs, identity)
    for i in range(m):
        assert sum(transform[i][k] * rhs[k] for k in range(m)) == b[i]
        for j, column in enumerate(columns):
            assert sum(transform[i][k] * column[k] for k in range(m)) == (d if i == j else 0)
    if m == n:
        solution = [F(v, d) for v in b]
        assert [sum(c[i] * x for c, x in zip(columns, solution)) for i in range(m)] == rhs


def test_residual_completion_preserves_original_equations():
    columns = [[1, 1, 0], [0, 1, 1]]
    addition = [0, 1, 0]
    rhs = [2, 5, 2]
    d, b, transformed = reduce_basis(columns, rhs, [addition])
    z = F(b[2], transformed[2][0])
    x = [(F(b[i])-transformed[i][0]*z)/d for i in range(2)]
    assert x == [2, 2] and z == 1
    assert [sum(c[i]*w for c,w in zip(columns,x))+addition[i]*z for i in range(3)] == rhs


def test_dependent_basis_refused():
    with pytest.raises(ValueError, match='Dependent'):
        reduce_basis([[1, 2], [2, 4]], [1, 2], [])


def test_noninteger_and_oversized_input_refused():
    with pytest.raises(ValueError, match='integer'):
        reduce_basis([[1.0]], [1], [])
    with pytest.raises(ValueError, match='integer'):
        reduce_basis([[1 << 4097]], [1], [])
