"""Exact codimension-one completion using supplied physical candidate columns."""
from fractions import Fraction as F
import numpy as np
from scipy.linalg import qr
from two_spectator_family_resume import solve_basis


def complete(matrix, rhs, candidates):
    m = len(rhs)
    if not 2 <= m <= 119 or len(matrix) != m or any(len(row) != m-1 for row in matrix):
        raise ValueError('One missing column in a bounded exact system required')
    if not 1 <= len(candidates) <= 5000 or any(len(col) != m for _, col in candidates):
        raise ValueError('Bounded complete physical candidate columns required')
    matrix = [list(map(F, row)) for row in matrix]
    rhs = list(map(F, rhs))
    # QR only nominates independent rows. The rational solve and null identity
    # below reject an incorrect nomination; floats never establish feasibility.
    _, _, pivots = qr(np.array(matrix, float).T, pivoting=True, mode='economic')
    keep = list(map(int, pivots[:-1]))
    omitted = int(pivots[-1])
    relation = solve_basis([[matrix[i][j] for i in keep] for j in range(m-1)],
                           [-matrix[omitted][j] for j in range(m-1)])
    left = [F(0)] * m
    left[omitted] = F(1)
    for i, value in zip(keep, relation):
        left[i] = value
    if any(sum(left[i]*matrix[i][j] for i in range(m)) for j in range(m-1)):
        raise ValueError('Exact left-null identity failed')
    residual = sum(a*b for a,b in zip(left,rhs))
    if not residual:
        raise ValueError('System already lies in the existing column span')
    possible = []
    for index, raw in candidates:
        col = list(map(F, raw))
        denominator = sum(a*b for a,b in zip(left,col))
        if denominator and residual/denominator > 0:
            possible.append((residual/denominator, index, col))
    possible.sort(key=lambda item:item[0])
    for attempt, (new_weight, index, col) in enumerate(possible[:4], 1):
        augmented = [row+[value] for row,value in zip(matrix,col)]
        weights = solve_basis(augmented,rhs)
        if weights[-1] != new_weight:
            raise ValueError('Exact completion identity failed')
        if min(weights) < 0:
            continue
        if any(sum(a*b for a,b in zip(row,weights)) != value for row,value in zip(augmented,rhs)):
            raise ValueError('Exact completed constraints failed')
        return index, weights, {'attempts':attempt,'eligible_physical_columns':len(possible),
                                'added_weight':str(new_weight),'exact_left_null_residual':str(residual)}
    raise ValueError('No nonnegative physical completion among the four bounded nominees')
