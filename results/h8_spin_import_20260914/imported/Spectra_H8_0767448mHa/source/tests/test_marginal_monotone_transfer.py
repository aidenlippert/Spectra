import copy
from fractions import Fraction as F
import unittest
from itertools import combinations
from unittest.mock import patch
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_monotone_transfer import replay
from experiments.marginal_symbolic import decode, encode, add, scale, mono, product
from experiments.marginal_spin_reduction import SpinZeroOracle, gap_coverage


def fixture(t=F(1, 6)):
    reference = build(4, 4, F(1, 3))
    target = build(4, 4, abs(t))
    if t < 0:
        h = decode(target['hamiltonian'], 8, 4)
        target['hamiltonian'] = encode({w: -v if len(w) == 2 else v for w, v in h.items()})
    return {'kind': 'joint_polynomial_monotone_transfer_v1', 'modes': 8,
            'particles': 4, 'hamiltonian': target['hamiltonian'],
            'target_lower': reference['target_lower'], 'reference_certificate': reference}


class MonotoneTransferTests(unittest.TestCase):
    def test_contraction_zero_and_sign_reversal_without_state_actions(self):
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action', side_effect=AssertionError('No states')):
            for t in [F(1, 3), F(1, 6), F(0), F(-1, 6), F(-1, 3)]:
                r = replay(fixture(t))
                self.assertEqual(r['determinant_actions'], 0)
                self.assertEqual(r['complement_lower'], '-1/1000')
                self.assertTrue(all(F(g['maximum_squared_magnitude_change']) <= 0 for g in r['changed_groups']))

    def test_target_may_have_a_sign_changing_spectator_amplitude(self):
        from experiments.marginal_polynomial_metric import JointPolynomial
        reference = build(4, 4, F(1, 3), slack=F(10))
        hopping = add(*(mono(((1, i), (0, j))) for i, j in [(0, 2), (2, 0), (1, 3), (3, 1)]))
        spectator = add(mono(((1, 4), (0, 4))), mono(((1, 5), (0, 5))))
        reference['hamiltonian'] = encode(add(decode(reference['hamiltonian'], 8, 4),
                                              scale(product(hopping, spectator), F(1, 12))))
        c = {'kind': 'joint_polynomial_monotone_transfer_v1', 'modes': 8, 'particles': 4,
             'hamiltonian': encode(add(decode(reference['hamiltonian'], 8, 4), scale(hopping, F(1, 4)))),
             'target_lower': reference['target_lower'], 'reference_certificate': reference}
        with self.assertRaisesRegex(ValueError, 'fixed-sign'):
            JointPolynomial(dict(c, polynomial_metric=reference['polynomial_metric'])).compile(F(c['target_lower']))
        result = replay(c)
        self.assertEqual(len(result['changed_groups']), 4)
        self.assertFalse(result['target_amplitude_sign_required'])

    def test_actual_target_reference_and_scope_refusals(self):
        with self.assertRaises(ValueError):
            replay(fixture(F(1, 2)))
        source = fixture()
        mutations = [lambda c: c.update(particles=2),
                     lambda c: c.update(target_lower='1'),
                     lambda c: c['reference_certificate'].update(kind='unverified_psd'),
                     lambda c: c['reference_certificate']['weight_proof'].update(bound=-1),
                     lambda c: c.update(hamiltonian=build(4, 3, F(1, 6))['hamiltonian'])]
        for mutation in mutations:
            c = copy.deepcopy(source)
            mutation(c)
            with self.assertRaises(ValueError):
                replay(c)
        c = fixture()
        h = decode(c['hamiltonian'], 8, 4)
        # A new edge has no reference amplitude and cannot be contractive.
        h = add(h, mono(((1, 0), (0, 6)), F(1, 100)), mono(((1, 6), (0, 0)), F(1, 100)),
                   mono(((1, 1), (0, 7)), F(1, 100)), mono(((1, 7), (0, 1)), F(1, 100)))
        c['hamiltonian'] = encode(h)
        with self.assertRaises(ValueError):
            replay(c)

    def test_energy_gap_dispatch_binds_outer_hamiltonian_and_full_valence(self):
        c = fixture()
        oracle = SpinZeroOracle(c)
        retained = [sum(1 << (2*i+(0 if i in alpha else 1)) for i in range(4))
                    for alpha in combinations(range(4), 2)]
        recipe = {'complement_atoms': {'kind': c['kind'], 'reference_certificate': c['reference_certificate']}}
        r = gap_coverage(oracle, retained, F(c['target_lower']), recipe)
        self.assertEqual(r['proof_family'], c['kind'])
        for hidden in ['hamiltonian', 'target_lower', 'modes']:
            bad = copy.deepcopy(recipe)
            bad['complement_atoms'][hidden] = c.get(hidden)
            with self.assertRaises(ValueError):
                gap_coverage(oracle, retained, F(c['target_lower']), bad)
        with self.assertRaises(ValueError):
            gap_coverage(oracle, retained[:-1], F(c['target_lower']), recipe)
        with self.assertRaises(ValueError):
            gap_coverage(SpinZeroOracle(build(4, 4, F(1, 2))), retained, F(c['target_lower']), recipe)


if __name__ == '__main__':
    unittest.main()
