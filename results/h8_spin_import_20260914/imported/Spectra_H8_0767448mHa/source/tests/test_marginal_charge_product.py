import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

from experiments.marginal_charge_product import feature_orbits, features, expand_recipe, replay, propose


def recipe():
    return {'kind': 'spin_charge_product_complement_v1', 'blocks': [{'states': [12]}],
            'metric_rule': {'family': 'pair_square', 'factors': ['1']*4, 'integer_scale': 10**8}}


class ChargeProductTests(unittest.TestCase):
    def test_pair_range_is_explicit_and_changes_feature_dictionary(self):
        self.assertEqual(len(feature_orbits(6, 'pair_square', 1)), 12)
        self.assertEqual(len(feature_orbits(6, 'pair_square', 2)), 16)
        self.assertEqual(feature_orbits(6, 'pair_square', 5), feature_orbits(6, 'pair_square'))
        for distance in (0, 6, True, 1.5):
            with self.assertRaises(ValueError): feature_orbits(6, 'pair_square', distance)
        c = recipe(); c['metric_rule']['max_pair_distance'] = 1
        self.assertEqual(expand_recipe(c, 4), expand_recipe(recipe(), 4))

    def test_small_exact_rounding_and_reflection_features(self):
        self.assertEqual(len(feature_orbits(6, 'pair')), 15)
        self.assertEqual(len(feature_orbits(6, 'pair_square')), 24)
        for s in range(16):
            reflected = ((s & 3) << 2) | (s >> 2)
            self.assertEqual(features(s, 2, 'pair_square'), features(reflected, 2, 'pair_square'))
        c = recipe(); c['metric_rule']['factors'] = ['2', '3', '5', '7']
        # q=(-1,+1): sum q=0, sum q²=2, q0*q1=-1, q0²*q1²=1.
        self.assertEqual(expand_recipe(c, 4)['blocks'][0]['metric_weights'], [1260000000])

    def test_direct_and_reference_energy_gap_are_rebuilt_and_bound(self):
        from tests.test_marginal_spin_temple import fixture
        from experiments.marginal_spin_temple import replay as replay_energy
        c = fixture(); red = c['spin_symmetric_certificate']; red.pop('blocks')
        red['complement_atoms'] = recipe()
        r = replay_energy(c)
        self.assertEqual(r['width'], replay_energy(fixture())['width'])
        self.assertEqual(r['complement_coverage']['proof_family'], recipe()['kind'])
        transferred = copy.deepcopy(c); red = transferred['spin_symmetric_certificate']
        atoms = red.pop('complement_atoms')
        red['complement_reference'] = {'hamiltonian': red['hamiltonian'],
                                      'complement_lower': red['complement_lower'], 'complement_atoms': atoms}
        self.assertEqual(replay_energy(transferred)['width'], r['width'])
        for mutate in (lambda red: red.update(complement_lower='5'),
                       lambda red: red['complement_atoms']['blocks'][0].update(states=[3])):
            bad = copy.deepcopy(c); mutate(bad['spin_symmetric_certificate'])
            with self.assertRaises(ValueError): replay_energy(bad)

    def test_malformed_metric_and_unphysical_coverage_rejected(self):
        for mutate in (lambda c: c['metric_rule'].update(factors=['1']),
                       lambda c: c['metric_rule']['factors'].__setitem__(0, '0'),
                       lambda c: c['metric_rule'].update(integer_scale=True),
                       lambda c: c['blocks'][0].update(metric_weights=[1]),
                       lambda c: c['blocks'][0].update(states=[True])):
            c = recipe(); mutate(c)
            with self.assertRaises(ValueError): expand_recipe(c, 4)

    def test_constructor_exports_replayable_compact_metric(self):
        from tests.test_marginal_spin_temple import fixture
        base = fixture()['spin_symmetric_certificate']
        base['blocks'] = [{'states': [12]}]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'source.json'; source.write_text(json.dumps(base))
            result = propose(source, root/'proof', target=F(7, 2))
            c = json.loads((root/'proof/certificate.json').read_text())
            self.assertEqual(replay(c), result)
            self.assertEqual(set(c['blocks'][0]), {'states'})
            self.assertEqual(result['metric_factor_count'], 4)
            with self.assertRaises(ValueError): propose(source, root/'proof', target=F(7, 2))


if __name__ == '__main__': unittest.main()
