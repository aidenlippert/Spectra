"""Exact identities, independent bit actions, and certificate refusal paths."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest

from experiments.marginal_symbolic import add, mono, product, scale, decode
from research.molecular_collective_20260913.core import density, matmul, transpose, eye
from research.certificate_scaling.commutator_dual_witness import psd
from research.response_consistency_20260913.separator import act
from research.compact_response_20260913 import program, closure, trial


def matrix_sum(*terms):
    return [[sum(c*a[i][j] for c, a in terms) for j in range(len(terms[0][1][0]))]
            for i in range(len(terms[0][1]))]


class Identities(unittest.TestCase):
    def test_compressed_density_square_includes_hopping(self):
        L = [[F(1), F(1, 3), F(1, 5)], [F(1, 3), F(-2, 7), F(-1, 4)],
             [F(1, 5), F(-1, 4), F(3, 2)]]
        operator = add(*(scale(density(i, j), L[i][j]) for i in range(3) for j in range(3)))
        core = add(mono((), 2*L[2][2]), *(scale(density(i, j), L[i][j]) for i in range(2) for j in range(2)))
        hopping = add(mono((), 2*sum(L[i][2]**2 for i in range(2))),
            *(scale(density(i, j), -L[i][2]*L[j][2]) for i in range(2) for j in range(2)))
        Q = program.projector(6)
        sandwich = lambda x: product(product(Q, x), Q)
        self.assertEqual(sandwich(product(operator, operator)), sandwich(add(product(core, core), hopping)))
        self.assertTrue(sandwich(hopping))

    def test_exact_noncommuting_chebyshev_response(self):
        D = [[F(3), F(1)], [F(1), F(2)]]; I = eye(2)
        B = [[F(1, 3), F(1, 5)], [F(1, 7), F(-1, 4)]]
        delta = F(1); M = F(4); c = (M+delta)/2; a = (M-delta)/2; z = c/a
        Z = matrix_sum((c/a, I), (-1/a, D)); k = 7
        Tprev, Tcur = I, Z
        qprev, qcur = matrix_sum((F(0), B)), matrix_sum((1/a, B))
        for j in range(1, k):
            qprev, qcur = qcur, matrix_sum((F(2), matmul(Z, qcur)), (F(-1), qprev), (2*program.chebyshev(j, z)/a, B))
            Tprev, Tcur = Tcur, matrix_sum((F(2), matmul(Z, Tcur)), (F(-1), Tprev))
        X = matrix_sum((1/program.chebyshev(k, z), qcur))
        E = matrix_sum((F(1), B), (F(-1), matmul(D, X)))
        exact_residual = matrix_sum((1/program.chebyshev(k, z), matmul(Tcur, B)))
        self.assertEqual(E, exact_residual)
        g = F(1); error_bound = matrix_sum((g*g/program.chebyshev(k, z)**2, I), (-1, matmul(transpose(E), E)))
        psd(error_bound)

    def test_rounded_response_congruence(self):
        A = [[12, 1], [1, 14]]; Bt = [[2, -1], [1, 3]]; D = [[8, 1], [1, 6]]
        den = 3; X = [[1, 2], [-1, 1]]; xden = 5; delta = F(1)
        K, kden, eta = closure.residual_and_K(A, Bt, D, den, X, xden, delta)
        H = [[F(x, den) for x in row] for row in [A[0]+[Bt[j][0] for j in range(2)],
            A[1]+[Bt[j][1] for j in range(2)], Bt[0]+D[0], Bt[1]+D[1]]]
        T = eye(4)
        for i in range(2):
            for j in range(2):
                T[2+i][j] = F(-X[i][j], xden)
        congruence = matmul(transpose(T), matmul(H, T))
        self.assertEqual([[F(x, kden) for x in row] for row in K], [row[:2] for row in congruence[:2]])
        self.assertEqual(eta, sum(congruence[i][j]**2 for i in range(2, 4) for j in range(2))/delta)

    def test_exact_factor_acceptance_and_failure(self):
        self.assertEqual(closure.check_factor([[5, 0], [0, 5]], 1, F(1), [[2], [0, 2]], 1), F(0))
        with self.assertRaises(ValueError):
            closure.check_factor([[4, 0], [0, 5]], 1, F(1), [[2], [0, 2]], 1)


class Certificates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.tail, cls.reference = program.load_case('h6')
        cls.response = json.loads((program.OUT/'h6_program.json').read_text())
        cls.cert = json.loads((program.OUT/'expanded_closure_certificate.json').read_text())
        cls.groups, cls.hden, cls.cost = closure.blocks(cls.data)

    def test_program_acceptance(self):
        for name in ('h6', 'h8'):
            data, tail, _ = program.load_case(name)
            cert = json.loads((program.OUT/f'{name}_program.json').read_text())
            receipt = program.check(data, tail, cert)
            self.assertLessEqual(F(receipt['exact_program_residual_penalty_Ha']), F(1, 1000000))
            self.assertEqual(receipt['many_body_states_enumerated'], 0)

    def test_program_false_bounds_and_binding_refused(self):
        for key, value in [('fixture_sha256', 'wrong'), ('tail_sha256', 'wrong'),
            ('delta_Ha', '1'), ('M_Ha', '10'), ('coupling_norm_Ha', '1'),
            ('order', 1), ('order', 257), ('order', True), ('delta_Ha', '0'),
            ('doubly_occupied_spatial_orbital', 0), ('target_Ha', 1.0)]:
            with self.subTest(key=key, value=value):
                cert = self.response | {key: value}
                with self.assertRaises(ValueError):
                    program.check(self.data, self.tail, cert)

    def test_complete_frozen_interval(self):
        receipt = closure.check(self.data, self.tail, self.response, self.cert, self.reference)
        self.assertEqual(F(receipt['width_Ha']), F(1, 1000))
        self.assertEqual(receipt['explicit_enumeration_and_matrices']['physical_sector_basis_labels_generated'], 924)

    def test_missing_and_duplicate_blocks_refused(self):
        for duplicate in (False, True):
            cert = deepcopy(self.cert)
            if duplicate:
                cert['blocks'][1] = cert['blocks'][0]
            else:
                cert['blocks'].pop()
            with self.assertRaises(ValueError):
                closure.check(self.data, self.tail, self.response, cert, self.reference)

    def test_response_tampering_refused(self):
        cert = deepcopy(self.cert)
        row = next(r for r in cert['blocks'] if r['response_columns'])
        row['response_columns'][0][0] += 10**15
        with self.assertRaisesRegex(ValueError, 'understated'):
            closure.check(self.data, self.tail, self.response, cert, self.reference)

    def test_factor_tampering_refused(self):
        cert = deepcopy(self.cert); cert['blocks'][0]['retained_factor'][0][0] += 10**14
        with self.assertRaises(ValueError):
            closure.check(self.data, self.tail, self.response, cert, self.reference)

    def test_closure_binding_and_denominator_refused(self):
        for key, value in [('response_program_sha256', 'wrong'), ('response_denominator', 0)]:
            with self.assertRaises(ValueError):
                closure.check(self.data, self.tail, self.response, self.cert | {key: value}, self.reference)

    def test_independent_car_actions_agree_with_blocks(self):
        h = decode(self.data['hamiltonian'], 12, 4)
        for group in self.groups:
            for j in sorted({0, len(group['states'])//2, len(group['states'])-1}):
                expected = act(h, {group['states'][j]: F(1)})
                actual = {s: F(group['H'][i][j], self.hden) for i, s in enumerate(group['states']) if group['H'][i][j]}
                self.assertEqual(actual, expected)

    def test_full_h8_expansion_refused(self):
        data, _, _ = program.load_case('h8')
        with self.assertRaisesRegex(ValueError, 'capped at H6'):
            closure.blocks(data)

    def test_scalar_obstruction(self):
        cert = json.loads((program.OUT/'scalar_obstruction.json').read_text())
        result = closure.check_scalar_obstruction(self.groups, self.hden, self.response, cert)
        self.assertLess(F(result['exact_expectation_Ha']), 0)
        with self.assertRaises(ValueError):
            closure.check_scalar_obstruction(self.groups, self.hden, self.response, cert | {'exact_expectation_Ha': '0'})

    def test_fresh_trial_acceptance_and_refusal(self):
        data, tail, reference = trial.load_case('fresh_h6_1p6')
        response = json.loads((program.OUT/'fresh_h6_1p6_program.json').read_text())
        cert = json.loads((program.OUT/'fresh_h6_1p6_trial.json').read_text())
        result = trial.check(data, tail, reference, response, cert)
        self.assertGreater(result['accepted_gain_mHa'], 0)
        with self.assertRaises(ValueError):
            trial.check(data, tail, reference, response, cert | {'upper_Ha': '0'})
        with self.assertRaises(ValueError):
            trial.check(data, tail, reference, response, cert | {'source_witness_sha256': 'wrong'})


if __name__ == '__main__':
    unittest.main()
