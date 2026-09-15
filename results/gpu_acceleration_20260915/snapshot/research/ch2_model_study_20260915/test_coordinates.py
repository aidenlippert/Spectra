import math
import unittest
from research.ch2_model_study_20260915.study import coordinates, chain_gradient, BOHR_ANGSTROM


class Coordinates(unittest.TestCase):
    def test_energy_gradient_chain_and_units(self):
        r, theta = 1.1, math.radians(113.)
        gradient = [[2., -1., .3], [.4, 2.3, -1.1], [.7, -1.6, .8]]

        def linear_energy(x):
            return sum(g*c/BOHR_ANGSTROM for (_, xyz), row in zip(coordinates(*x), gradient) for g, c in zip(row, xyz))

        analytic = chain_gradient(r, theta, gradient)
        for j in range(2):
            plus, minus = [r, theta], [r, theta]
            plus[j] += 1e-6
            minus[j] -= 1e-6
            self.assertAlmostEqual(analytic[j], (linear_energy(plus)-linear_energy(minus))/2e-6, places=8)


if __name__ == '__main__':
    unittest.main()
