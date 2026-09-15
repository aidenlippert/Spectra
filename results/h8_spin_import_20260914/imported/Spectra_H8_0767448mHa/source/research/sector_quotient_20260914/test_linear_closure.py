from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import number_shift, product, mono, scale, add, canonical, adj


class LinearClosureTests(unittest.TestCase):
    def test_right_number_ideal_for_annihilator_and_adjoint(self):
        for m, n in ((6, 3), (12, 6), (16, 8)):
            for i in range(m):
                a = mono(((0, i),)); creation = canonical(adj(a))
                D = product(number_shift(m, 0), a)
                self.assertEqual(add(D, scale(a, -(n-1)), scale(product(a, number_shift(m, n)), -1)), {})
                self.assertEqual(add(canonical(adj(D)), scale(creation, -n), scale(product(creation, number_shift(m, n)), -1)), {})


if __name__ == '__main__': unittest.main()
