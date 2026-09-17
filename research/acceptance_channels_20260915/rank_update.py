"""Exact-form low-rank normal update; floating point is proposal-only."""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import LinearOperator


class RankUpdatedNormal(LinearOperator):
    def __init__(self,base,columns):
        self.base=base;self.columns=columns
        super().__init__(dtype=np.dtype('float64'),shape=base.shape)

    def _matvec(self,value):
        return self.base@value+self.columns@(self.columns.T@value)

    def __add__(self,other):
        if sparse.issparse(other):return RankUpdatedNormal(self.base+other,self.columns)
        return super().__add__(other)


def factorize(matrix,splu):
    factor=splu(matrix.base.tocsc())
    W=factor.solve(matrix.columns.toarray())
    small=np.eye(matrix.columns.shape[1])+matrix.columns.T@W
    class Factor:
        def solve(self,rhs):
            initial=factor.solve(rhs)
            return initial-W@np.linalg.solve(small,matrix.columns.T@initial)
    return Factor()
