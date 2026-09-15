"""Differential checks of the projection contract, including failure inputs."""
import os
import unittest
import numpy as np
from research.gpu_acceleration_20260915.kernel_bench import cpu_project, GPUProject, HybridProject


class ProjectionTests(unittest.TestCase):
    def matrices(self):
        rng = np.random.default_rng(17)
        result = [np.zeros((4,4)), np.diag([-2., 0., 0., 3.]), -np.eye(12), np.eye(6)]
        for n in (2, 6, 36, 90, 222, 260):
            a = rng.normal(size=(n,n))
            result.append((a+a.T)/2)
        return result

    def compare(self, project):
        matrices = self.matrices()
        expected = cpu_project(matrices)
        actual = project(matrices)
        self.assertEqual(len(actual), len(expected))
        for q, got, want in zip(matrices, actual, expected):
            np.testing.assert_allclose(got, want, atol=2e-10, rtol=2e-10)
            np.testing.assert_allclose(got, got.T, atol=2e-10)
            self.assertGreaterEqual(np.linalg.eigvalsh(got)[0], -2e-10)
            self.assertGreaterEqual(np.linalg.eigvalsh(got-q)[0], -2e-10)
        with self.assertRaises(ValueError):
            project([np.full((4,4), np.nan)])

    def test_cpu_evd(self):
        self.compare(lambda blocks: cpu_project(blocks, 'evd'))

    def test_cpu_positive(self):
        self.compare(lambda blocks: cpu_project(blocks, 'positive'))

    @unittest.skipUnless(os.environ.get('SPECTRA_TEST_GPU') == '1', 'requires explicit CUDA test')
    def test_gpu_mixed(self):
        gpu = GPUProject([q.shape for q in self.matrices()], mixed=True)
        self.compare(gpu)

    @unittest.skipUnless(os.environ.get('SPECTRA_TEST_GPU') == '1', 'requires explicit CUDA test')
    def test_hybrid(self):
        self.compare(HybridProject([q.shape for q in self.matrices()]))


if __name__ == '__main__':
    unittest.main()
