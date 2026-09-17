"""Independent small matrix and analytic oracles for the accepting boundary."""
from fractions import Fraction as F
import copy
import math
import unittest

from experiments.marginal_symbolic import encode, mono, add
from research.intervention_reduction_20260916.exact import (
    check, comparison_envelope, digest, exact_moments, integrate_envelope, inverse_spd,
)


def matrix_poly(a):
    return add(*(mono(((1,i),(0,j)), F(x)) for i,row in enumerate(a) for j,x in enumerate(row) if x))


def case(h, vectors, controls=()):
    data = {'modes':len(h), 'particles':1, 'hamiltonian':encode(matrix_poly(h))}
    proposal = {'kind':'integer_control_subspace_v1', 'fixture_sha256':digest(data),
                'configurations':[1<<i for i in range(len(vectors))], 'vectors':vectors,
                'denominator':1,
                'controls':[{'operator':encode(matrix_poly(c)), 'amplitude_Ha':str(a)} for c,a in controls]}
    return data, proposal


class ExactReductionTest(unittest.TestCase):
    def test_nonorthogonal_projection_and_external_action(self):
        h = [[2,1,3], [1,4,2], [3,2,7]]
        data, p = case(h, [[1,0],[1,1]])
        m = exact_moments(data,p)
        self.assertEqual(m['metric'], [[F(2),F(1)],[F(1),F(1)]])
        # Direct elementary vector arithmetic: V0=(1,1,0), V1=(0,1,0).
        self.assertEqual(m['projected'][0], [[F(8),F(5)],[F(5),F(4)]])
        self.assertEqual(m['generators'][0], [[F(3),F(1)],[F(2),F(3)]])
        self.assertEqual(m['residual_diagonal_blocks'][0][0], [F(25),F(4)])
        self.assertEqual(m['counts']['external_configurations'],1)

    def test_control_cross_terms_cannot_be_discarded(self):
        h = [[0,1], [1,0]]
        data,p = case(h, [[1]], [(h,F(1))])
        m = exact_moments(data,p)
        _, radii, squares = comparison_envelope(m)
        self.assertEqual(squares,[F(4)])
        self.assertEqual(radii,[F(2)])
        r,_ = check(data,p,horizon=F(1,10))
        self.assertGreaterEqual(F(r['state_vector_error_bound']),F(1,5))

    def test_diagonal_phase_has_zero_transition_envelope(self):
        h = [[10000,0],[0,-9000]]
        data,p = case(h, [[1,0],[0,1]])
        r,_ = check(data,p)
        self.assertEqual(r['comparison_matrix'],[['0','0'],['0','0']])
        self.assertEqual(F(r['state_vector_error_bound']),0)

    def test_envelope_taylor_matches_analytic_integral(self):
        c = [[F(0),F(1,4)],[F(1,4),F(0)]]
        bound,tail = integrate_envelope(c,[F(0),F(3)],F(2),order=8)
        exact_numeric = 12*(math.cosh(.5)-1)
        self.assertGreaterEqual(float(bound)+1e-14,exact_numeric)
        self.assertLess(float(bound)-exact_numeric,1e-8)
        self.assertGreater(tail,0)

    def test_singular_metric_and_wrong_bindings_refuse(self):
        data,p = case([[1,0],[0,2]], [[1,2],[1,2]])
        with self.assertRaisesRegex(ValueError,'positive definite'):
            exact_moments(data,p)
        p['fixture_sha256'] = 'wrong'
        with self.assertRaisesRegex(ValueError,'binding'):
            exact_moments(data,p)

    def test_asymmetric_and_wrong_sector_refuse(self):
        with self.assertRaisesRegex(ValueError,'Hermitian'):
            data,p = case([[0,1],[0,0]], [[1]])
            exact_moments(data,p)
        data,p = case([[1,0],[0,2]], [[1]])
        p['configurations']=[3]
        with self.assertRaisesRegex(ValueError,'fixed-N'):
            exact_moments(data,p)

    def test_inverse_independent_multiplication(self):
        a=[[F(3),F(1),F(1)],[F(1),F(4),F(2)],[F(1),F(2),F(5)]]
        b=inverse_spd(a)
        self.assertEqual([[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)],
                         [[F(i==j) for j in range(3)] for i in range(3)])

    def test_invalid_rational_control_and_tail_order_refuse(self):
        data,p=case([[1,0],[0,2]],[[1]], [([[0,1],[1,0]],F(1))])
        bad=copy.deepcopy(p);bad['controls'][0]['amplitude_Ha']=0.1
        with self.assertRaisesRegex(ValueError,'strings'):
            exact_moments(data,bad)
        with self.assertRaisesRegex(ValueError,'order'):
            integrate_envelope([[F(10)]],[F(1)],F(10),0)


if __name__ == '__main__':
    unittest.main()
