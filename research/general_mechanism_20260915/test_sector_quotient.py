"""Independent small-Fock-space oracles and rejection tests for the quotient."""
from fractions import Fraction as F
from itertools import combinations
import unittest

from experiments.marginal_symbolic import add, mono, product
from research.general_mechanism_20260915.sector_quotient import (
    operator_quotient, psd_quotient,
)


def act(operator, state):
    result = {}
    for bits, amplitude in state.items():
        for word, coefficient in operator.items():
            value, occupation = coefficient*amplitude, bits
            for creation, mode in reversed(word):
                occupied = (occupation >> mode) & 1
                if occupied == creation:
                    value = 0
                    break
                value *= (-1)**((occupation & ((1 << mode)-1)).bit_count())
                occupation ^= 1 << mode
            result[occupation] = result.get(occupation, 0)+value
    return {k: v for k, v in result.items() if v}


def singlets_two_electrons(spatial):
    states = [{(1 << (2*p)) | (1 << (2*p+1)): F(1)} for p in range(spatial)]
    for p, q in combinations(range(spatial), 2):
        states.append({(1 << (2*p)) | (1 << (2*q+1)): F(1),
                       (1 << (2*p+1)) | (1 << (2*q)): F(-1)})
    return states


def dictionary(spatial):
    number = add(*(mono(((1, i), (0, i))) for i in range(2*spatial)))
    raising = add(*(mono(((1, 2*p), (0, 2*p+1))) for p in range(spatial)))
    a = mono(((0, 1),))
    return [a, product(a, add(number, mono((), -2))), product(raising, a),
            mono(((0, 3),)), mono(((1, 0),)), mono(((1, 2),)),
            mono(()), raising, {}, product(mono(((1, 2), (0, 2))), a)]


class SectorQuotientTest(unittest.TestCase):
    def test_complete_physical_action_quotient_against_independent_oracle(self):
        for spatial in (2, 3):
            operators = dictionary(spatial)
            states = singlets_two_electrons(spatial)
            result = operator_quotient(operators, 2*spatial, 2)
            self.assertEqual(result['sector_dimension'], len(states))
            self.assertEqual(result['determinants_enumerated'], 0)
            by_index = {i: (g, k) for g in result['groups']
                        for k, i in enumerate(g['indices'])}
            # Compare every Gram entry, including all omitted cross-group ones.
            for i, left in enumerate(operators):
                for j, right in enumerate(operators):
                    exact = F(0)
                    for state in states:
                        x, y = act(left, state), act(right, state)
                        exact += sum(v*y.get(k, 0) for k, v in x.items())/sum(v*v for v in state.values())
                    exact /= len(states)
                    group, a = by_index[i]
                    other, b = by_index[j]
                    got = group['gram'][a][b] if group is other else 0
                    self.assertEqual(got, exact)
            for group in result['groups']:
                local = [operators[i] for i in group['indices']]
                for vector in group['kernel']:
                    for state in states:
                        combined = {}
                        for coefficient, operator in zip(vector, local):
                            for bits, value in act(operator, state).items():
                                combined[bits] = combined.get(bits, 0)+coefficient*value
                        self.assertFalse({k: v for k, v in combined.items() if v})
            # Singlet-null is not zero as an operator on the entire Fock space.
            self.assertTrue(act(operators[2], {(1 << 1) | (1 << 3): F(1)}))
            self.assertFalse(result['existing_SOS_ideal_equivalence_proved'])

    def test_particle_and_spin_conditions_are_not_ignored(self):
        operators = dictionary(3)
        self.assertEqual(operator_quotient([operators[1]], 6, 2)['rank'], 0)
        self.assertEqual(operator_quotient([operators[1]], 6, 4)['rank'], 1)
        with self.assertRaises(ValueError):
            operator_quotient(operators, 6, 3)

    def test_reject_inexact_or_unsupported_dictionary(self):
        cases = [{(): 0.1}, {((0, 4),): F(1)},
                 {((0, 0),): F(1), ((0, 2),): F(1)},
                 {((2, 0),): F(1)}, {((1, 0),)*4: F(1)}]
        for case in cases:
            with self.assertRaises(ValueError):
                operator_quotient([case], 4, 2)
        with self.assertRaises(ValueError):
            operator_quotient([], 4, 2)

    def test_exact_psd_and_kernel_gates(self):
        result = psd_quotient([[1, 2, 0], [2, 4, 0], [0, 0, 0]])
        self.assertEqual(result['pivots'], [0])
        self.assertEqual(result['kernel'], [[F(-2), F(1), F(0)], [F(0), F(0), F(1)]])
        for bad in ([[1, 2], [2, 1]], [[0, 1], [1, 0]], [[1, 1], [0, 1]], [[1.0]], [[1, 2]]):
            with self.assertRaises(ValueError):
                psd_quotient(bad)


if __name__ == '__main__':
    unittest.main()
