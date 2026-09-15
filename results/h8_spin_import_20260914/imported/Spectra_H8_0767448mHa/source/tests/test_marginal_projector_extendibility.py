import json
import unittest
from fractions import Fraction as F
from pathlib import Path
from math import lcm

from experiments.marginal_projector_extendibility import projector_bound, replay
from experiments.marginal_polynomial_sos import integer_psd


class ProjectorExtendibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (Path(__file__).resolve().parents[1]
                / 'results/marginal_graded_hubbard8/weighted_window_family/certificate.json')
        cls.vector = json.loads(path.read_text())['upper_vector']

    def _independent_two_window_columns(self):
        amplitudes = {int(s): a for s, a in self.vector.items() if a}
        norm = sum(a * a for a in amplitudes.values())
        # Direct tensor products in base-four site occupation order, without
        # the production routine's general offset/environment decomposition.
        columns = [{s+256*exterior:a for s,a in amplitudes.items()}
                   for exterior in range(4)]
        columns += [{exterior+4*s:a for s,a in amplitudes.items()}
                    for exterior in range(4)]
        return columns, norm

    def test_two_window_exact_ceiling_and_rejects_too_small(self):
        result = projector_bound(self.vector, 2, '1.380638468')
        self.assertTrue(result['accepted'])
        self.assertEqual(F(result['average_fidelity_ceiling']), F(690319234, 10**9))
        self.assertEqual(result['gram_dimension'], 8)
        self.assertEqual(result['maximum_psd_dimension'], 2)
        columns, norm = self._independent_two_window_columns()
        gram = [[sum(columns[i].get(s, 0) * columns[j].get(s, 0)
                     for s in set(columns[i]) | set(columns[j]))
                 for j in range(8)] for i in range(8)]
        # Independent sparse tensor embedding: the cross overlap is diagonal.
        overlap = [[F(gram[b][4 + a], norm) for a in range(4)] for b in range(4)]
        self.assertEqual(overlap, [[F(106854722802607664801035,
                                      895219092079179430700602), 0, 0, 0],
                                    [0, F(170377411618491025274633,
                                          447609546039589715350301), 0, 0],
                                    [0, 0, F(170377411618491025274633,
                                             447609546039589715350301), 0],
                                    [0, 0, 0, F(106854722802607664801035,
                                                895219092079179430700602)]])
        ceiling = F('1.380638468')
        rational = [[(ceiling * norm if i == j else 0) - gram[i][j] for j in range(8)]
                    for i in range(8)]
        scale = lcm(*(x.denominator for row in rational for x in row))
        self.assertTrue(integer_psd([[int(x*scale) for x in row] for row in rational])
                        ['positive_semidefinite'])
        with self.assertRaises(ValueError):
            projector_bound(self.vector, 2, '1.38')

    def test_three_window_bound_and_product_state_edge(self):
        result = projector_bound(self.vector, 3, '1.901418483')
        self.assertEqual(result['gram_dimension'], 48)
        self.assertEqual(F(result['average_fidelity_ceiling']), F(1901418483, 3 * 10**9))
        product = projector_bound({85: 1}, 2, 2)
        self.assertEqual(product['average_fidelity_ceiling'], '1')
        self.assertTrue(product['spin_sector_split'])
        result4 = projector_bound(self.vector, 4, '2.223786409')
        self.assertEqual(result4['gram_dimension'], 256)
        with self.assertRaises(ValueError):
            projector_bound(self.vector, 4, '2.22')

    def test_replay_valid_penalized_certificate(self):
        certificate = {
            'kind': 'hubbard_projector_extension_v1',
            'a': '531373/1000000', 'b': '3/4',
            'vector': self.vector, 'windows': 2,
            'projector_sum_ceiling': '1.380638468',
            'penalty': '0.3157', 'penalized_lower': '-1.724725',
            'chain_sites': 8,
        }
        result = replay(certificate)
        self.assertTrue(result['accepted'])
        self.assertEqual(result['local_sum_dimensions'], 256)
        self.assertLess(F(result['periodic_lower_density']), F('-0.647'))
        self.assertGreater(F(result['periodic_lower_density']), F('-2.040424674') / 3)
        self.assertEqual(F(result['open_lower_density']),
                         F(result['periodic_lower_density']) - F(2, 8))
        self.assertEqual(F(result['open_lower_energy']),
                         F(result['open_lower_density']) * 8)

    def test_replay_refuses_unsound_or_malformed_certificates(self):
        base = {
            'kind': 'hubbard_projector_extension_v1',
            'a': '531373/1000000', 'b': '3/4',
            'vector': self.vector, 'windows': 2,
            'projector_sum_ceiling': '1.380638468',
            'penalty': '0.3157', 'penalized_lower': '-1.724725',
            'chain_sites': 8,
        }
        for field, value in [('projector_sum_ceiling', '1.38'),
                             ('penalty', '-1/10'),
                             ('penalized_lower', '-1.7'),
                             ('chain_sites', 6), ('windows', 1)]:
            with self.subTest(field=field):
                bad = dict(base, **{field: value})
                with self.assertRaises(ValueError):
                    replay(bad)
        with self.assertRaises(ValueError):
            replay(dict(base, vector={85: 1}))
        with self.assertRaises(ValueError):
            replay(dict(base, vector={'15': 1.0}))

    def test_reflection_and_spin_projector_gates(self):
        base = {
            'kind': 'hubbard_projector_extension_v1', 'a': '531373/1000000',
            'b': '3/4', 'vector': self.vector, 'windows': 2,
            'projector_sum_ceiling': '1.380638468', 'penalty': '0.3157',
            'penalized_lower': '-1.724725', 'chain_sites': 8,
        }
        with self.assertRaisesRegex(ValueError, 'spin-number'):
            replay(dict(base, vector={15: 1, 85: 1}))
        # A single four-particle basis state is not reflection parity eigenstate.
        with self.assertRaisesRegex(ValueError, 'reflection'):
            replay(dict(base, vector={15: 1}))


if __name__ == '__main__':
    unittest.main()
