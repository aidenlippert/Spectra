import unittest
import numpy as np

from .mps_numeric import inner, norm, linear_combination, apply_mpo, compress


def product_mps(vecs):
    return [np.asarray(v, dtype=float).reshape(1, len(v), 1) for v in vecs]


def dense(s):
    x = s[0]
    for a in s[1:]: x = np.einsum("...i,ijk->...jk", x, a).reshape(-1, a.shape[2])
    return x[:, 0]


class MPSNumericTests(unittest.TestCase):
    def test_inner_and_norm(self):
        a = product_mps(([1, 2], [3, -1], [2, 1]))
        b = product_mps(([2, 0], [1, 4], [1, -2]))
        self.assertAlmostEqual(inner(a, b), float(np.vdot(dense(a), dense(b))))
        self.assertAlmostEqual(norm(a), np.linalg.norm(dense(a)))

    def test_direct_sum_linear_combination(self):
        a = product_mps(([1, 2], [3, 1]))
        b = product_mps(([2, -1], [1, 4]))
        c = linear_combination([a, b], [2.0, -0.5])
        self.assertTrue(np.allclose(dense(c), 2 * dense(a) - .5 * dense(b)))

    def test_apply_product_mpo(self):
        a = product_mps(([1, 2], [3, 1]))
        # local X on site zero, diagonal Z on site one
        X = np.array([[0., 1.], [1., 0.]])
        Z = np.diag([1., -1.])
        W = [X.reshape(1, 1, 2, 2), Z.reshape(1, 1, 2, 2)]
        got = apply_mpo(a, W)
        want = np.kron(X, Z) @ dense(a)
        self.assertTrue(np.allclose(dense(got), want))

    def test_compression_reports_error_and_matches_dense(self):
        rng = np.random.default_rng(4)
        s = [rng.normal(size=(1, 2, 3)), rng.normal(size=(3, 2, 4)),
             rng.normal(size=(4, 2, 1))]
        old = dense(s)
        c, d = compress(s, 2)
        self.assertLessEqual(c[0].shape[2], 2)
        self.assertGreaterEqual(d["discarded_squared_norm"], 0.)
        self.assertAlmostEqual(d["discarded_squared_norm"], np.linalg.norm(old) ** 2 - np.linalg.norm(dense(c)) ** 2, places=8)


if __name__ == "__main__": unittest.main()
