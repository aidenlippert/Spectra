from fractions import Fraction as F
from itertools import product
import unittest
from .tensor_operator import build, element, product_terms, expand_products


def literal_action(word, ket):
    sign = 1
    for creation, mode in reversed(word):
        occupied = bool(ket & (1 << mode))
        if occupied == bool(creation):
            return None, 0
        sign *= (-1) ** sum(bool(ket & (1 << j)) for j in range(mode))
        ket ^= 1 << mode
    return ket, sign


class TestOperatorMPO(unittest.TestCase):
    def test_short_words_independent_oracle(self):
        words = [()] + [w for n in (1, 2, 3) for w in product(tuple(product((0, 1), range(2))), repeat=n)]
        words += [((1, 1), (1, 0), (0, 1), (0, 0)), ((1, 0), (1, 0))]
        for word in words:
            data = {'modes': 2, 'hamiltonian': [{'word': word, 'coefficient': '7/13'}]}
            op = build(data)
            self.assertEqual(expand_products(op), product_terms(data))
            for ket in range(4):
                out, sign = literal_action(word, ket)
                for bra in range(4):
                    self.assertEqual(element(op, bra, ket), F(7 * sign, 13) if out == bra else F(0))

    def test_mixed_quartic_and_exact_reconstruction(self):
        words = [((1, 3), (1, 0), (0, 2), (0, 1)),
                 ((1, 1), (1, 2), (0, 0), (0, 3)),
                 ((1, 0), (0, 0)), ((0, 2), (1, 2)), ()]
        coeffs = [F(7, 11), F(-5, 17), F(2, 3), F(-3, 5), F(19, 23)]
        data = {'modes': 4, 'hamiltonian': [
            {'word': w, 'coefficient': str(c)} for w, c in zip(words, coeffs)]}
        op = build(data)
        self.assertEqual(expand_products(op), product_terms(data))
        for ket in range(16):
            expected = {}
            for word, c in zip(words, coeffs):
                out, sign = literal_action(word, ket)
                if sign:
                    expected[out] = expected.get(out, 0) + sign * c
            for bra in range(16):
                self.assertEqual(element(op, bra, ket), expected.get(bra, F(0)))

    def test_invalid_schema_refused(self):
        for word in [[(2, 0)], [(1, 4)], [('1', 0)]]:
            with self.assertRaises((ValueError, TypeError)):
                build({'modes': 4, 'hamiltonian': [{'word': word, 'coefficient': '1'}]})


if __name__ == '__main__':
    unittest.main()
