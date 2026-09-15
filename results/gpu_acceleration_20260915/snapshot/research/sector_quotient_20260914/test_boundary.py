"""Small independent analytic cone cases validate the numerical proposer."""
import unittest
import numpy as np
from scipy import sparse
from research.sector_quotient_20260914.boundary import solve


class BoundaryTests(unittest.TestCase):
    def test_smallest_eigenvalue_with_free_bound(self):
        H = np.array([[2., .4], [.4, -1.]])
        def A(Q):
            return Q[0].ravel()
        def AT(y):
            Y = y.reshape(2, 2)
            return [(Y+Y.T)/2]
        free = sparse.csc_matrix(np.eye(2).ravel()[:, None])
        Q, x, y, Z, stats = solve(A, AT, free, H.ravel(), [H+2*np.eye(2)], np.array([-2.]), 5,
                                  max_outer=1000, tolerance=1e-9)
        self.assertAlmostEqual(x[0], np.linalg.eigvalsh(H)[0], places=7)
        self.assertLess(np.linalg.norm(A(Q)+free@x-H.ravel()), 1e-7)

    def test_two_scalar_cones_and_an_unrestricted_multiplier(self):
        # h=(2,0,5); equations Q1+b+t=2, t=0, Q2+2b=5.
        # Both Q>=0, hence max b=2.
        def A(Q):
            return np.array([Q[0][0, 0], 0., Q[1][0, 0]])
        def AT(y):
            return [np.array([[y[0]]]), np.array([[y[2]]])]
        free = sparse.csc_matrix([[1., 1.], [0., 1.], [2., 0.]])
        Q, x, y, Z, stats = solve(A, AT, free, np.array([2., 0., 5.]),
                                  [np.array([[2.]]), np.array([[5.]])], np.zeros(2), 5,
                                  max_outer=1000, tolerance=1e-9)
        self.assertAlmostEqual(x[0], 2., places=7)
        self.assertAlmostEqual(x[1], 0., places=7)


if __name__ == '__main__':
    unittest.main()
