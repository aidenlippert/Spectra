"""Independent redundant-coordinate checks against the existing dense projector."""
import unittest
import numpy as np
from scipy import sparse
from research.sector_quotient_20260914.eliminated import Quotient
from research.nvidia_followup_20260915.sparse_quotient import SparseQuotient


class Problem:
    def __init__(self):
        rng = np.random.default_rng(835)
        base = rng.normal(size=(31, 5))
        self.free = sparse.csc_matrix(np.column_stack([rng.normal(size=31),base,base[:,0]+2*base[:,3],base[:,4],np.zeros(31)]))
        self.rhs = rng.normal(size=31)
        self.M = rng.normal(size=(31, 9))

    def A(self,Q):
        return self.M@Q[0].ravel()

    def AT(self,y):
        return [(self.M.T@y).reshape(3,3)]


class SparseTests(unittest.TestCase):
    def test_same_projector_and_recovered_coefficients(self):
        op=Problem()
        a,b=Quotient(op),SparseQuotient(op)
        self.assertEqual(a.record['numerical_ideal_rank'],b.record['numerical_ideal_rank'])
        rng=np.random.default_rng(826)
        for _ in range(4):
            y=rng.normal(size=31)
            np.testing.assert_allclose(a.project(y),b.project(y),atol=2e-12)
            Q=[rng.normal(size=(3,3))]
            ax,bx=a.recover(Q),b.recover(Q)
            self.assertAlmostEqual(ax[0],bx[0],places=12)
            np.testing.assert_allclose(op.free@ax,op.free@bx,atol=2e-12)

    def test_dependent_energy_coordinate_refused(self):
        op=Problem()
        data=op.free.toarray()
        data[:,0]=data[:,1]
        op.free=sparse.csc_matrix(data)
        with self.assertRaises(ValueError):
            SparseQuotient(op)


if __name__=='__main__':
    unittest.main()
