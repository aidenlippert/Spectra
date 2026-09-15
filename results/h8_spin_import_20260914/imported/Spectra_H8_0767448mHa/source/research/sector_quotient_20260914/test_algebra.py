"""Exact tests; the 64-input bit oracle is an algebra diagnostic only."""
from fractions import Fraction as F
from itertools import combinations
import unittest
from research.sector_quotient_20260914.algebra import (
    lift, contract, inverse_contracted, combine, split_three_body, inner,
    factorized_lift_check, from_polynomial, to_polynomial)
from experiments.marginal_symbolic import product, number_shift, add, scale


def action(poly, bits):
    out = {}
    for word, coefficient in poly.items():
        state = bits
        for creation, orbital in reversed(word):
            occupied = (state >> orbital) & 1
            if occupied == creation:
                coefficient = 0
                break
            coefficient *= (-1)**((state & ((1 << orbital)-1)).bit_count())
            state ^= 1 << orbital
        if coefficient:
            out[state] = out.get(state, 0) + coefficient
    return {key: value for key, value in out.items() if value}


class QuotientTests(unittest.TestCase):
    def test_all_325_matrix_units(self):
        count = 0
        for m in (5, 6):
            pairs = list(combinations(range(m), 2))
            for I in pairs:
                for J in pairs:
                    V = {(I, J): F(1)}
                    W = lift(V, m, 2)
                    K = contract(W, m, 3)
                    self.assertEqual(K, combine((m-4, V), (1, lift(contract(V, m, 2), m, 1))))
                    self.assertEqual(inverse_contracted(K, m), V)
                    self.assertEqual(split_three_body(W, m), (V, {}))
                    count += 1
        self.assertEqual(count, 325)

    def test_sparse_12_and_16(self):
        for m in (12, 16):
            W = {((0, 2, 3), (1, 2, m-1)): F(7, 9), ((0, 1, 2), (0, 1, 2)): F(-5, 11)}
            V, Z = split_three_body(W, m)
            self.assertEqual(contract(Z, m, 3), {})
            self.assertEqual(inner(lift(V, m, 2), Z), 0)
            self.assertEqual(split_three_body(Z, m), ({}, Z))

    def test_contraction_zero_is_not_zero(self):
        W = {((0, 1, 2), (3, 4, 5)): F(1), ((3, 4, 5), (0, 1, 2)): F(1)}
        self.assertEqual(contract(W, 6, 3), {})
        self.assertEqual(split_three_body(W, 6), ({}, W))
        factors = [(F(1, 2), {(0, 1, 2): 1, (3, 4, 5): 1}),
                   (F(-1, 2), {(0, 1, 2): 1, (3, 4, 5): -1})]
        V, receipt = factorized_lift_check(factors, 6)
        self.assertEqual(V, {})
        self.assertEqual(receipt['remainder_squared_norm'], 2)
        self.assertFalse(receipt['is_lift'])

    def test_nonzero_positive_completion(self):
        V = {(I, I): F(-1) if I == (0, 1) else F(1, 2) for I in combinations(range(6), 2)}
        W = lift(V, 6, 2)
        factors = [(v, {I: 1}) for (I, J), v in W.items()]
        recovered, receipt = factorized_lift_check(factors, 6)
        self.assertEqual(recovered, V)
        self.assertTrue(receipt['is_lift'])
        self.assertGreater(receipt['squared_norm'], 0)

    def test_nonzero_completion_through_unchanged_exact_acceptor(self):
        from experiments.marginal_symbolic import encode, verified_residual
        V = {(I, I): F(-1) if I == (0, 1) else F(1, 2) for I in combinations(range(6), 2)}
        W = lift(V, 6, 2)
        self.assertTrue(all(v == F(3, 2) for v in W.values()))
        # 3/2 = 1^2 + (1/2)^2 + (1/2)^2, so no irrational factors.
        blocks = [{'words': [[(0, i) for i in reversed(I)]], 'factor': [[2], [1], [1]]}
                  for (I, J), value in W.items()]
        p = to_polynomial(V, 6, 2)
        core = {'modes': 6, 'particles': 3, 'operator_degree': 3,
                'hamiltonian': encode(p), 'b': '0', 'number_multiplier': encode(scale(p, -1)),
                'denominator': 2, 'blocks': blocks}
        residual, receipt = verified_residual(core)
        self.assertEqual(residual, {})
        self.assertEqual(F(receipt['lower']), 0)
        core['particles'] = 2
        wrong_residual, wrong_receipt = verified_residual(core)
        self.assertTrue(wrong_residual)
        self.assertLessEqual(F(wrong_receipt['lower']), -1)

    def test_independent_bit_action_and_multiplier_sign(self):
        m = 6
        V = {((0, 1), (2, 4)): F(3, 7), ((2, 4), (0, 1)): F(3, 7), ((0, 3), (0, 3)): F(-2, 5)}
        p = to_polynomial(V, m, 2)
        w = to_polynomial(lift(V, m, 2), m, 3)
        self.assertEqual(from_polynomial(p, m, 2), V)
        self.assertEqual(product(number_shift(m, 2), p), w)
        for bits in range(1 << m):
            number = bits.bit_count()
            self.assertEqual(action(w, bits), action(scale(p, number-2), bits))
            # H=(N-2)V, S=L(V), X=-V: H-S-(Nhat-N)X=0.
            residual = add(scale(p, number-2), scale(w, -1), product(number_shift(m, number), p))
            self.assertEqual(residual, {})

    def test_refusal(self):
        for bad in ({((0, 0), (0, 1)): F(1)}, {((0, 1), (0, 6)): F(1)}, {((0, 1), (0, 1)): .5}):
            with self.assertRaises(ValueError):
                lift(bad, 6, 2)
        with self.assertRaises(ValueError):
            split_three_body({}, 4)
        with self.assertRaises(ValueError):
            factorized_lift_check([(1., {(0, 1, 2): 1})], 6)


if __name__ == '__main__':
    unittest.main()
