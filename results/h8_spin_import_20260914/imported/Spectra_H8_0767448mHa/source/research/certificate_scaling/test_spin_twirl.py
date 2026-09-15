"""Exact Schur-formula controls for the SU(2) polynomial twirl."""
import unittest
from fractions import Fraction as F

from experiments.marginal_hunt_car import adj
from experiments.marginal_symbolic import add, product, scale,number_shift
from research.certificate_scaling.spin_basis import spin_action
from research.certificate_scaling.spin_twirl import twirl
from math import comb, lcm
from experiments.marginal_spin_reduction import spin_operators,commutator


def descend(p, j2):
    out = [p]
    for r in range(1, j2 + 1):
        out.append(scale(spin_action(out[-1], raising=False), F(1, r)))
    return out


def square(p):
    return product(adj(p), p)


class SpinTwirlControls(unittest.TestCase):
    def check_formula(self, j2, copies, coeffs):
        descendants = [descend(p, j2) for p in copies]
        wden = lcm(*(comb(j2, r) for r in range(j2 + 1)))
        weights = [wden // comb(j2, r) for r in range(j2 + 1)]
        p = {}
        for a, ds in enumerate(descendants):
            for r, q in enumerate(ds):
                p = add(p, scale(q, coeffs[a][r]))
        actual = twirl(square(p))
        expected = {}
        for a in range(len(copies)):
            for b in range(len(copies)):
                for r in range(j2 + 1):
                    factor = coeffs[a][r] * coeffs[b][r] / F((j2 + 1) * weights[r])
                    invariant = {}
                    for s in range(j2 + 1):
                        invariant = add(invariant, scale(product(adj(descendants[a][s]), descendants[b][s]), weights[s]))
                    expected = add(expected, scale(invariant, factor))
        self.assertEqual(actual, expected)
        self.assertEqual(twirl(actual), actual)
        modes=2*((max(i for w in p for _,i in w)//2)+1)
        plus,minus,z=spin_operators(modes)
        for generator in (plus,minus,z):self.assertFalse(commutator(generator,actual))

    def test_nonorthogonal_multiplicity_copies_half_and_one(self):
        self.check_formula(1, [
            {((1, 0),): F(1)},
            {((1, 0),): F(1), ((1, 2),): F(2)},
        ], [[F(2), F(3)], [F(-1), F(2)]])
        self.check_formula(2, [
            {((1, 0), (1, 2)): F(1)},
            {((1, 0), (1, 2)): F(1), ((1, 0), (1, 4)): F(2)},
        ], [[F(1), F(2), F(-1)], [F(2), F(-1), F(3)]])

    def test_three_halves_and_singlet(self):
        self.check_formula(3, [{((1, 0), (1, 2), (1, 4)): F(1)}], [[F(1), F(-2), F(3), F(1)]])
        singlet = {((1, 0), (1, 1)): F(1)}
        self.assertEqual(twirl(square(singlet)), square(singlet))

    def test_rejects_odd_or_too_high_degree(self):
        with self.assertRaises(ValueError): twirl({((1, 0),): F(1)})
        with self.assertRaises(ValueError): twirl({tuple((1, i) for i in range(7)): F(1)})

    def test_wrong_descendant_weights_change_the_operator(self):
        p={((1,0),(1,2)):F(1)}
        wrong=scale(add(*(square(q) for q in descend(p,2))),F(1,3))
        self.assertNotEqual(twirl(square(p)),wrong)

    def test_number_ideal_commutes_with_twirl(self):
        x={((1,0),(0,0)):F(2,3),((1,0),(1,2),(0,0),(0,2)):F(-1,7)}
        shift=number_shift(4,2)
        self.assertEqual(twirl(product(shift,x)),product(shift,twirl(x)))


if __name__ == '__main__':
    unittest.main()
