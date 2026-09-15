"""Exact acceptance, response algebra and metric regression checks."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest

from research.compact_response_20260913 import program
from research.global_response_20260913 import global_program as g
from research.global_response_20260913.reference_diagnostic import scalar_terminal_margin


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def add(a, b, scale=F(1)):
    return [[x+scale*y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def mul(c, a):
    return [[c*x for x in row] for row in a]


class GlobalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.tail, _ = g.load_case('h6')
        cls.cert = json.loads((g.OUT/'h6_global.json').read_text())

    def test_valid_global_certificate_and_partition(self):
        r = g.check(self.data, self.tail, self.cert)
        self.assertEqual(r['partition'], 'union_last_two_double')
        self.assertEqual(r['H_actions_per_K_vector'], 237)
        self.assertFalse(r['terminal_positivity_proved_here'])
        self.assertEqual(r['many_body_states_enumerated'], 0)

    def test_reject_changed_input_and_target_type(self):
        c = deepcopy(self.cert)
        c['fixture_sha256'] = '0'*64
        with self.assertRaises(ValueError): g.check(self.data, self.tail, c)
        c = deepcopy(self.cert)
        c['target_Ha'] = float(F(c['target_Ha']))
        with self.assertRaises(ValueError): g.check(self.data, self.tail, c)

    def test_reject_forged_response_and_partition(self):
        c = deepcopy(self.cert)
        c['response']['order'] -= 1
        with self.assertRaises(ValueError): g.check(self.data, self.tail, c)
        c = deepcopy(self.cert)
        c['first_sector']['orbital'] -= 1
        with self.assertRaises(ValueError): g.check(self.data, self.tail, c)

    def test_nonpositive_denominator_and_allowance_refused(self):
        bound = {'delta': F(-1), 'M': F(5), 'g': F(1)}
        self.assertEqual(g.scalar_program(bound, F(1, 10**6))['status'], 'no_certified_denominator')
        with self.assertRaises(ValueError): g.scalar_program(bound, F(0))

    def test_matrix_chebyshev_factorization(self):
        # A nondiagonal rational D tests the actual operator identity, rather
        # than a commuting replacement of the surrounding P/H/Q products.
        eye = [[F(1), F(0)], [F(0), F(1)]]
        d = [[F(3), F(1, 3)], [F(1, 3), F(4)]]
        inv = mul(F(1)/(d[0][0]*d[1][1]-d[0][1]**2),
                  [[d[1][1], -d[0][1]], [-d[1][0], d[0][0]]])
        z = mul(F(1, 3), add(mul(F(4), eye), d, F(-1)))
        def p(k):
            prev, cur = eye, z
            for _ in range(1, k): prev, cur = cur, add(mul(2, mm(z, cur)), prev, F(-1))
            return mm(add(eye, mul(1/program.chebyshev(k, F(4, 3)), cur), F(-1)), inv)
        for k in (1, 2, 4):
            pk = p(k)
            left = add(mul(2, pk), mm(mm(pk, d), pk), F(-1))
            t = program.chebyshev(2*k, F(4, 3))
            self.assertEqual(left, mul(t/(t+1), p(2*k)))

    def test_nonunitary_error_requires_metric(self):
        # T=diag(1/10,1), M=diag(-1,0): T*MT >= -.01 I,
        # yet M is not >= -.01 I. The inverse norm pays exactly 100.
        t = [[F(1, 10), F(0)], [F(0), F(1)]]
        m = [[F(-1), F(0)], [F(0), F(0)]]
        transformed = mm(mm(t, m), t)
        epsilon = F(1, 100)
        self.assertEqual(transformed[0][0], -epsilon)
        self.assertLess(m[0][0], -epsilon)
        self.assertEqual(m[0][0], -epsilon*100)
        relative = add(transformed, mul(epsilon, mm(t, t)))
        self.assertLess(relative[0][0], 0)  # relative certificate correctly fails

    def test_terminal_allowances_paid_once(self):
        r = scalar_terminal_margin(F(1), F(0), [F(1, 10), F(1, 5)])
        self.assertEqual(F(r['terminal_margin_Ha']), F(7, 10))
        self.assertFalse(scalar_terminal_margin(F(1), F(1), [])['positive'])
        self.assertFalse(scalar_terminal_margin(F(1), F(0), [F(2)])['positive'])
        with self.assertRaises(ValueError): scalar_terminal_margin(F(1), F(0), [F(-1)])


if __name__ == '__main__': unittest.main()
