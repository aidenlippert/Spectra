import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_molecular_gershgorin import SOURCES, matrix, upper, schur_pivots, replay
from experiments.marginal_symbolic import decode, encode, add, mono

ROOT = Path(__file__).resolve().parents[1]


class MolecularGershgorinTests(unittest.TestCase):
    def test_matrix_matches_independent_ordered_occupation_oracle(self):
        for name, path in SOURCES.items():
            c = json.loads((ROOT / path).read_text())
            states, a = matrix(c)
            self.assertEqual(states, sorted(states))
            oracle = [[F(0)] * 70 for _ in range(70)]
            for word, coefficient in decode(c['hamiltonian'], 8, 4).items():
                for j, state in enumerate(states):
                    occupied = [i for i in range(8) if state >> i & 1]
                    sign = 1
                    for creation, mode in reversed(word):
                        if bool(creation) == (mode in occupied):
                            break
                        position = sum(i < mode for i in occupied)
                        sign *= (-1) ** position
                        if creation: occupied.insert(position, mode)
                        else: occupied.pop(position)
                    else:
                        destination = sum(1 << i for i in occupied)
                        oracle[states.index(destination)][j] += sign * coefficient
            self.assertEqual(a, oracle, name)
            expected = F(c['independent_upper']['upper'])
            self.assertEqual(upper(a, c['independent_upper']), expected)
            changed = copy.deepcopy(c['independent_upper'])
            changed['upper'] = '1000'
            self.assertEqual(upper(a, changed), expected)

    def test_negative_schur_denominator_cannot_pass_positive_retained_matrix(self):
        data = ([[F(2)]], [[F(1)]], F(0))
        # Omitting the denominator gate would accept b=1 incorrectly.
        self.assertIsNone(schur_pivots(data, F(1)))
        self.assertIsNone(schur_pivots(data, F(0)))
        self.assertIsNotNone(schur_pivots(data, F(-1)))

    def test_saved_certificates_and_false_witnesses(self):
        for name in ('rectangle', 'square', 'rectangle_coupling', 'square_coupling',
                     'rectangle_witness', 'square_witness'):
            p = ROOT / 'results/marginal_molecular_gershgorin' / name
            c = json.loads((p / 'certificate.json').read_text())
            got = replay(c)
            expected = json.loads((p / 'receipt.json').read_text())
            self.assertTrue(all(got[k] == expected[k] for k in ('lower', 'upper', 'width')))
            self.assertGreater(F(got['positive_denominator']), 0)
        for field, value in (('retained_states', [15, 15]), ('retained_states', [0]),
                             ('retained_states', [True]), ('modes', 10), ('particles', 3),
                             ('lower', got['complement_lower']),
                             ('independent_upper', {'amplitudes': [0] * 70})):
            bad = copy.deepcopy(c)
            bad[field] = value
            with self.assertRaises(ValueError): replay(bad)
        bad = copy.deepcopy(c)
        bad['hamiltonian'] = encode(add(decode(c['hamiltonian'], 8, 4), mono((), F(-1))))
        with self.assertRaises(ValueError): replay(bad)


if __name__ == '__main__':
    unittest.main()
