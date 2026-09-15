"""Numerical coordinate elimination checks, independent of molecular results."""
from types import SimpleNamespace
import unittest
import numpy as np
from scipy import sparse
from research.sector_quotient_20260914.eliminated import Quotient


class EliminationTests(unittest.TestCase):
    def test_duplicate_ideal_coordinates_and_bound_recovery(self):
        F = sparse.csc_matrix([[1., 1., 2.], [0., 1., 2.], [2., 0., 0.]])
        op = SimpleNamespace(free=F, rhs=np.array([3., 1., 4.]),
            A=lambda Q: np.array([Q[0][0, 0], 0., Q[1][0, 0]]),
            AT=lambda y: [np.array([[y[0]]]), np.array([[y[2]]])])
        q = Quotient(op)
        self.assertEqual(q.record['numerical_ideal_rank'], 1)
        self.assertTrue(np.allclose(np.column_stack([q.project(v) for v in F.toarray().T]), 0, atol=1e-12))
        self.assertTrue(np.allclose(F.T@q.y0, [1., 0., 0.], atol=1e-12))
        Q = [np.array([[1.]]), np.array([[2.]])]
        recovered = q.recover(Q)
        self.assertAlmostEqual(recovered[0], 1., places=12)
        self.assertTrue(np.allclose(F@recovered+op.A(Q), op.rhs, atol=1e-12))


if __name__ == '__main__':
    unittest.main()
