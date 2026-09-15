import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_charge_product import feature_orbits, features
from experiments.marginal_coherent_tree import CoherentCharge, replay, build
from experiments.marginal_spin_constructor import spin_states
from experiments.marginal_symbolic import mono, add, product, adj, scale, encode


def fixture():
    n = [mono(((1, i), (0, i))) for i in range(8)]
    h = add(*(scale(product(n[i], n[i+1]), 4) for i in range(0, 8, 2)))
    hop = add(*(mono(((1, spin), (0, 2+spin)), F(-1, 3)) for spin in (0, 1)))
    hop = add(hop, adj(hop))
    h = add(h, hop, scale(product(hop, add(n[4], n[5])), F(-2, 3)))
    pair = mono(((1, 0), (1, 1), (0, 3), (0, 2)), F(1, 5))
    h = add(h, pair, adj(pair))
    count = len(feature_orbits(4, 'pair_square', 1))
    return {'modes': 8, 'particles': 4, 'hamiltonian': encode(h),
            'metric_rule': {'family': 'pair_square', 'max_pair_distance': 1,
                            'factors': [str(F(10+i, 11+i)) for i in range(count)]}}


def direct_rows(o):
    states = spin_states(o.oracle)
    qstates = [s for s in states if any(((s >> (2*i)) & 3) in (0, 3) for i in range(o.sites))]
    factors = [factor for _, factor in o.local_factors]
    def weight(s):
        q = [((s >> (2*i)) & 3).bit_count()-1 for i in range(o.sites)]
        out = F(1)
        for label, factor in o.local_factors:
            exponent = 1
            for i, power in label: exponent *= q[i]**power
            out *= factor**exponent
        return out
    weights = {s: weight(s) for s in qstates}
    rows = {}
    for s in qstates:
        row = o.oracle.action(s)
        rows[s] = row.get(s, F(0))-sum((abs(v)*weights[t]/weights[s]
                    for t, v in row.items() if t != s and t in weights), F(0))
    return rows


class CoherentTreeTests(unittest.TestCase):
    def test_tree_gap_binds_complete_valence_reference_in_energy_and_transfer(self):
        import contextlib, io
        from tests.test_marginal_valence_reference import fixture as small_fixture
        from experiments.marginal_valence_reference import build as build_energy
        from experiments.marginal_spin_temple import replay as replay_energy
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'source.json'; source.write_text(json.dumps(small_fixture()))
            with contextlib.redirect_stdout(io.StringIO()): build_energy(source, root/'proof', upper_steps=2)
            c = json.loads((root/'proof/certificate.json').read_text()); expected = replay_energy(c)
            red = c['spin_symmetric_certificate']; rule = {'family': 'pair_square', 'factors': ['1']*4}
            o = CoherentCharge(dict(red, metric_rule=rule)); tree, _ = o.cover(F(red['complement_lower']))
            recipe = {'kind': 'valence_coherent_tree_v1', 'metric_rule': rule, 'tree': tree, 'max_nodes': 10000}
            red['complement_atoms'] = recipe
            self.assertEqual(replay_energy(c)['width'], expected['width'])
            self.assertEqual(replay_energy(c)['complement_coverage']['determinant_actions'], 0)
            transferred = copy.deepcopy(c); r = transferred['spin_symmetric_certificate']; atoms = r.pop('complement_atoms')
            r['complement_reference'] = {'hamiltonian': r['hamiltonian'], 'complement_lower': r['complement_lower'], 'complement_atoms': atoms}
            self.assertEqual(replay_energy(transferred)['width'], expected['width'])
            for mutate in (lambda r: r.update(retained_states=[r['retained_states'][0]]),
                           lambda r: r.update(retained_states=[3, 12]),
                           lambda r: r.update(complement_lower='100'),
                           lambda r: r['complement_atoms'].update(hamiltonian=r['hamiltonian'])):
                bad = copy.deepcopy(c); mutate(bad['spin_symmetric_certificate'])
                with self.assertRaises(ValueError): replay_energy(bad)

    def test_source_indicator_polynomial_preserves_exclusivity_and_partial_bounds(self):
        data = fixture(); data.update(source_polynomial_bound=True, diagonal_charge_lambda='1')
        o = CoherentCharge(data); rows = direct_rows(o)
        for depth in (0, 2, 4, 6, 8):
            mask = (1 << depth)-1
            for bits in {s & mask for s in rows}:
                closed = o.close(mask, bits)
                bound = o.branch_lower(*closed)
                self.assertLessEqual(bound, min(v for s, v in rows.items() if s & mask == bits))
        # One alpha particle occupies one of two available modes. A negative
        # coefficient on their product must contribute zero, not a penalty.
        mask = o.all_bits ^ ((1 << 0) | (1 << 2))
        bits = (1 << 4) | (1 << 1) | (1 << 3)
        polynomial = {1 << 0: F(2), 1 << 2: F(3), (1 << 0) | (1 << 2): F(-100)}
        self.assertEqual(o.cardinality_polynomial_lower(polynomial, mask, bits), 2)
        tree, r = o.cover(F(2))
        self.assertEqual(o.cover(F(2), tree)[1]['covered_Q_configurations'], 30)
        self.assertGreater(r['maximum_source_polynomial_terms'], 0)

    def test_charge_quadratic_tangents_bound_diagonal_and_reject_false_threshold(self):
        data = fixture(); data['diagonal_charge_lambda'] = '1'
        o = CoherentCharge(data); rows = direct_rows(o)
        self.assertGreater(o.charge_diagonal_lower(0, 0), o.oracle.diagonal_lower(0, 0))
        for depth in (0, 2, 4, 6, 8):
            mask = (1 << depth)-1
            for bits in {s & mask for s in rows}:
                closed = o.close(mask, bits)
                physical = [o.oracle.action(s).get(s, F(0)) for s in rows if s & mask == bits]
                self.assertLessEqual(o.charge_diagonal_lower(*closed), min(physical))
                self.assertLessEqual(o.branch_lower(*closed), min(v for s, v in rows.items() if s & mask == bits))
        data['diagonal_charge_lambda'] = '100'
        with self.assertRaises(ValueError): CoherentCharge(data)

    def test_coherent_rows_equal_physical_rows_with_density_and_pair_hopping(self):
        o = CoherentCharge(fixture()); rows = direct_rows(o)
        self.assertEqual(len(rows), 30)
        for s, lower in rows.items():
            self.assertEqual(o.branch_lower(o.all_bits, s), lower)

    def test_partial_occupation_bounds_are_conservative(self):
        o = CoherentCharge(fixture()); rows = direct_rows(o)
        for depth in (0, 2, 4, 6):
            mask = (1 << depth)-1
            for bits in {s & mask for s in rows}:
                closed = o.close(mask, bits)
                self.assertIsNotNone(closed)
                actual = min(value for s, value in rows.items() if s & mask == bits)
                self.assertLessEqual(o.branch_lower(*closed), actual)

    def test_build_replay_no_actions_and_invalid_coverage_refusals(self):
        data = fixture()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'source.json'; source.write_text(json.dumps(data))
            with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                       side_effect=AssertionError('No determinant actions')):
                receipt = build(source, root/'proof', F(2))
                c = json.loads((root/'proof/certificate.json').read_text())
                self.assertEqual(replay(c), receipt)
            self.assertEqual(receipt['covered_Q_configurations'], 30)
            self.assertEqual(receipt['determinant_actions'], 0)
            self.assertLessEqual(receipt['exact_occupation_leaves'], 30)
            for update in ({'tree': None}, {'tree': ['valence']}, {'tree': ['split', True, [], []]},
                           {'target_lower': '100', 'tree': ['bound']}, {'max_nodes': 0}):
                bad = copy.deepcopy(c); bad.update(update)
                with self.assertRaises(ValueError): replay(bad)
            with self.assertRaises(ValueError): build(source, root/'proof', F(2))

    def test_combinatorial_sector_coverage_without_enumeration_at_64_modes(self):
        from math import comb
        data = {'kind': 'valence_coherent_tree_v1', 'modes': 64, 'particles': 32,
                'hamiltonian': encode(add(*(mono(((1, i), (0, i))) for i in range(64)))),
                'metric_rule': {'family': 'pair_square', 'max_pair_distance': 1,
                                'factors': ['1']*len(feature_orbits(32, 'pair_square', 1))},
                'target_lower': '32', 'tree': ['bound'], 'max_nodes': 1}
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                   side_effect=AssertionError('No determinant actions')):
            r = replay(data)
        self.assertEqual(r['covered_Q_configurations'], comb(32, 16)**2-comb(32, 16))
        self.assertEqual(r['exact_occupation_leaves'], 0)
        self.assertEqual(r['nodes'], 1)


if __name__ == '__main__': unittest.main()
