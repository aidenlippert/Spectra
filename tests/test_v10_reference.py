import unittest
from fractions import Fraction as F
from itertools import product
from copy import deepcopy

import numpy as np

from experiments.v10_reference import (
    BudgetExceeded, Reference, at_zero, check_expansion, construct, hermitian_modes,
    pauli_to_eigen,eigen_to_pauli,evaluate_expansion,
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
            for word, c in ref.v.items():
                q = np.array([[1]], complex)
                for letter in word:
                    q = np.kron(q, mats[letter])
                he += float(c) * q
            for word in map(''.join,product('IZPM',repeat=n)):
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
        layers, _ = construct(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1, 500), max_order=3)
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
        h = {'Z': F(1),'X':F(1)}
        layers, _ = construct(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1, 100), max_order=3)
        missing=deepcopy(layers)
        # Add a real Hermitian constant: differential and initial tests must
        # reject it even though conjugate-pair validation succeeds.
        missing[1][('Z',F(0),F(0),0)]=(F(1),F(0))
        self.assertTrue(hermitian_modes(missing[1]))
        with self.assertRaises(ValueError): check_expansion(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1), missing)
        bad = dict(layers[0]); key = next(iter(bad)); bad[key] = (bad[key][0], bad[key][1] + 1)
        self.assertFalse(hermitian_modes(bad))
        with self.assertRaises(ValueError): check_expansion(h, F(1, 2), 1, {'X': F(1)}, F(1, 5), F(1), [bad])
        wrong=deepcopy(layers)
        wrong[1]={k:(2*c[0],2*c[1]) for k,c in wrong[1].items()}
        self.assertTrue(hermitian_modes(wrong[1]));self.assertEqual(at_zero(wrong[1]),{})
        with self.assertRaises(ValueError):check_expansion(h,F(1,2),1,{'X':F(1)},F(1,5),F(1),wrong)

    def test_basis_round_trip(self):
        op={p:F(i-20,17) for i,p in enumerate(map(''.join,product('IXYZ',repeat=3)))}
        back=eigen_to_pauli(pauli_to_eigen(op,3),3)
        self.assertEqual(back,{p:(c,F(0)) for p,c in op.items() if c})

    def test_complete_two_qubit_output_against_dense_unitary(self):
        # An independent physical endpoint: diagonalize the 4x4 Hermitian H,
        # rather than differentiating the same symbolic recurrence again.
        I=np.eye(2);X=np.array([[0,1],[1,0]],complex)
        Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1])
        mats={'I':I,'X':X,'Y':Y,'Z':Z}
        h={'ZI':F(1,2),'IZ':F(-1,3),'XX':F(2,5),'ZZ':F(1,7)}
        o={'ZI':F(1)};T=F(1,5)
        layers,_=construct(h,0,2,o,T,F(1,10**7))
        receipt=check_expansion(h,0,2,o,T,F(1,10**7),layers)
        out,err,_=evaluate_expansion(layers,2,T,F(1,10**9))
        self.assertEqual(receipt['status'],'certified')
        dense=lambda p:np.kron(mats[p[0]],mats[p[1]])
        H=sum(float(c)*dense(p) for p,c in h.items())
        vals,vecs=np.linalg.eigh(H)
        U=(vecs*np.exp(1j*float(T)*vals))@vecs.conj().T
        expected=U@dense('ZI')@U.conj().T
        actual=sum(float(c)*dense(p) for p,c in out.items())
        delta=np.linalg.norm(expected-actual,2)
        self.assertLessEqual(delta,float(F(receipt['bound'])+err)+1e-13)
        self.assertGreater(delta,1e-14)

    def test_rejects_caps(self):
        with self.assertRaises(BudgetExceeded):
            Reference({'X': F(1)}, F(0), 1, cap=1).column('Z')
        with self.assertRaises(BudgetExceeded):
            construct({'Z': F(1)}, F(1), 1, {'X': F(1)}, F(1), F(1, 100), max_order=2, cap=1)


if __name__ == '__main__':
    unittest.main()
