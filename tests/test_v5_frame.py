import unittest
from fractions import Fraction as F
from itertools import product
from experiments.v5_frame import Q, BASE_LAMBDAS, frame, null_probe, covariance, dot


class FrameTests(unittest.TestCase):
    def test_dense_orthogonal_frame_all_hypotheses(self):
        self.assertTrue(all(v != 0 for col in Q for v in col))
        for i, a in enumerate(Q):
            for j, b in enumerate(Q): self.assertEqual(dot(a, b), int(i == j))
        for signs in product((-1, 1), repeat=3):
            modes, r = frame(signs)
            vectors = modes + (r,)
            for i, a in enumerate(vectors):
                self.assertTrue(all(v != 0 for v in a))
                for j, b in enumerate(vectors): self.assertEqual(dot(a, b), int(i == j))

    def test_null_probe_all_prefixes_and_current_alternatives(self):
        for length in range(3):
            for prefix in product((-1, 1), repeat=length):
                probe = null_probe(prefix)
                self.assertEqual(dot(probe, probe), 1)
                modes, _ = frame(prefix)
                self.assertTrue(all(dot(probe, mode) == 0 for mode in modes))
                plus, _ = frame(prefix + (1,))
                minus, _ = frame(prefix + (-1,))
                self.assertEqual(dot(probe, plus[-1]), 0)
                self.assertEqual(abs(dot(probe, minus[-1])), F(24, 25))

    def test_dense_covariance_eigenrelations(self):
        for signs in product((-1, 1), repeat=3):
            modes, residual = frame(signs)
            matrix = covariance(signs, F(3, 2))
            self.assertTrue(all(matrix[i][j] != 0 for i in range(4) for j in range(4)))
            for vector, eigenvalue in zip(modes + (residual,), tuple(1 + F(3, 2) * l for l in BASE_LAMBDAS) + (1,)):
                actual = tuple(sum(matrix[i][j] * vector[j] for j in range(4)) for i in range(4))
                self.assertEqual(actual, tuple(eigenvalue * x for x in vector))

    def test_validation(self):
        for value in ((True,), (0,), (1, 1, 1, 1), [1]):
            with self.assertRaises(ValueError): frame(value)
        with self.assertRaises(ValueError): covariance((1,), 1)
        with self.assertRaises(ValueError): covariance((), F(1))
        with self.assertRaises(ValueError): null_probe((1, 1, 1))


if __name__ == '__main__': unittest.main()
