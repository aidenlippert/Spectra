import unittest
from fractions import Fraction as F

import numpy as np

from experiments.v10_reference import (
    BudgetExceeded, Reference, at_zero, check_expansion, construct, hermitian_modes,
)


class V10ReferenceTests(unittest.TestCase):
    def test_local_column_matches_dense_commutator_n1_and_n2(self):
        # Compare every eigenword column against direct 2^n x 2^n matrices.
        I = np.eye(2, dtype=complex)
        X = np.array([[0, 1], [1, 0]], complex)
        Y = np.array([[0, -1j], [1j, 0]], complex)
        Z = np.diag([1, -1]).astype(complex)
        mats = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
        for n, h in ((1, {'Z': F(2)}), (2, {'ZI': F(1), 'IX': F(2), 'XY': F(1, 2)})):
            ref = Reference(h, F(1, 3), n)
            he = np.zeros((2 ** n, 2 ** n), complex)
            for word, c in h.items():
                q = np.array([[1]], complex)
                for letter in word:
                    q = np.kron(q, mats[letter])
                he += float(c) * q
            for word in ref.ve:
                got = ref.column(word)
                actual = np.zeros_like(he)
                q = np.array([[1]], complex)
                for letter in word:
                    q = np.kron(q, {'I': I, 'Z': Z, 'P': (X + 1j * Y) / 2,
                                    'M': (X - 1j * Y) / 2}[letter])
                expected = 1j * (he @ q - q @ he)
                for outword, coeff in got.items():
                    r = np.array([[1]], complex)
                    for letter in outword:
                        r = np.kron(r, {'I': I, 'Z': Z, 'P': (X + 1j * Y) / 2,
                                        'M': (X - 1j * Y) / 2}[letter])
                    actual += complex(float(coeff[0]), float(coeff[1])) * r
                np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_construct_and_complete_recurrence(self):
        h = {'Z': F(1), 'X': F(1)}
        layers, _ = construct(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1, 10**12), max_order=3)
        result = check_expansion(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1), layers)
        self.assertEqual(result['status'], 'certified')
        self.assertEqual(result['order'], 3)

    def test_repeated_frequency_and_zero_initial_condition(self):
        ref = Reference({'X': F(1)}, F(0), 1)
        forcing = ref.apply_b(ref.initial({'Z': F(1)}))
        integrated = ref.integrate(forcing)
        self.assertTrue(any(k == 1 for _, _, _, k in integrated))
        self.assertEqual(at_zero(integrated), {})

    def test_rejects_missing_initial_hermiticity_and_recurrence_tampering(self):
        h = {'Z': F(1)}
        layers, _ = construct(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1, 10**12), max_order=2)
        missing = list(layers)
        missing[1] = dict(missing[1]); key = next(iter(missing[1])); missing[1][key] = (missing[1][key][0] + 1, missing[1][key][1])
        with self.assertRaises(ValueError): check_expansion(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1), missing)
        bad = dict(layers[0]); key = next(iter(bad)); bad[key] = (bad[key][0], bad[key][1] + 1)
        self.assertFalse(hermitian_modes(bad))
        with self.assertRaises(ValueError): check_expansion(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1), [bad])

    def test_rejects_caps(self):
        with self.assertRaises(BudgetExceeded):
            Reference({'X': F(1)}, F(0), 1, cap=1).column('Z')
        with self.assertRaises(BudgetExceeded):
            construct({'Z': F(1)}, F(1), 1, {'X': F(1)}, F(1), F(1, 100), max_order=2, cap=1)


if __name__ == '__main__':
    unittest.main()
