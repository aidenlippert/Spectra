import unittest
import numpy as np
from scipy import sparse
from research.collective_completion_20260914.moment_dual import problem


class MomentDualTests(unittest.TestCase):
    def test_dual_squares_reconstruct_with_offdiagonals_and_scaling(self):
        # H - b I = Q, with the deliberately nonunit scaling used in search.
        P=sparse.eye(3,format='csc');free=sparse.csc_matrix([[1.],[0.],[1.]])
        rhs=np.array([2.,2.,3.]);ii,jj=np.triu_indices(2)
        p,eq,y,cones=problem(free,rhs,[(P,ii,jj)],np.array([.3,2.,.7]))
        p.solve(solver='CLARABEL')
        x=-eq.dual_value;Q=cones[0].dual_value
        np.testing.assert_allclose(free@x+P@Q[ii,jj],rhs,atol=1e-7)
        self.assertGreater(np.linalg.eigvalsh(Q)[0],-1e-8)
        self.assertAlmostEqual(float(x[0]),float(np.linalg.eigvalsh([[2,2],[2,3]])[0]),places=7)


if __name__=='__main__':unittest.main()
