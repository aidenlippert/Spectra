import unittest
import numpy as np
from research.acceptance_channels_20260915.fixed_scalar import step


class FixedScalarTest(unittest.TestCase):
    def test_affine_psd_intersection_on_known_matrix(self):
        # Solve Q + z I = H-bI, with non-diagonal positive target.
        H = np.array([[2., .3], [.3, 1.]])
        F = np.eye(2).reshape(4, 1)
        G = np.eye(4)+F@F.T
        target = (H-.5*np.eye(2)).ravel()
        A = lambda blocks: blocks[0].ravel()
        AT = lambda y: [y.reshape(2, 2)]
        def project(blocks):
            ev, V = np.linalg.eigh((blocks[0]+blocks[0].T)/2)
            return [(V*np.maximum(0., ev))@V.T]
        W, z = [np.array([[0., 1.], [1., -2.]])], np.zeros(1)
        for _ in range(100):
            W, z = step(A, AT, F, target, W, z, lambda r: np.linalg.solve(G, r), project)
        Q = project(W)
        np.testing.assert_allclose(A(Q)+F@z, target, atol=1e-10)
        self.assertGreaterEqual(np.linalg.eigvalsh(Q[0]).min(), -1e-12)


if __name__ == '__main__':
    unittest.main()
