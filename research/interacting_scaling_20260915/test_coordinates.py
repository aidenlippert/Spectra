from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from experiments.marginal_symbolic import add, canonical, decode, encode, mono, number_shift
from research.collective_completion_20260914.spin_screen import spin_squared
from research.interacting_scaling_20260915.rotation import rotate_polynomial
from research.interacting_scaling_20260915.coupling import coupled_fixture


class CoordinateTests(unittest.TestCase):
    Z = [[3, 0, 4, 0], [0, 5, 0, 0], [-4, 0, 3, 0], [0, 0, 0, 5]]

    def test_charge_spin_and_full_degree_family_transport(self):
        for poly in [number_shift(8, 4), spin_squared(8)]:
            self.assertEqual(rotate_polynomial(poly, self.Z, 5), canonical(poly))
        p = mono(((1, 0), (0, 3), (0, 2)))
        q = rotate_polynomial(p, self.Z, 5)
        self.assertTrue(all(len(w) == 3 and sum(2*c-1 for c, i in w) == -1 for w in q))
        transpose = [list(row) for row in zip(*self.Z)]
        self.assertEqual(rotate_polynomial(q, transpose, 5), canonical(p))

    def test_old_index_restricted_triples_are_not_rotation_invariant(self):
        p = mono(((0, 0), (0, 1), (0, 2)))
        q = rotate_polynomial(p, self.Z, 5)
        self.assertTrue(any(any(i < 4 for c, i in w) and any(i >= 4 for c, i in w) for w in q))

    def test_complete_coupling_endpoint_and_charge_transfer(self):
        h = add(mono(((1, 0), (0, 0))), mono(((1, 0), (0, 4)), F(-1, 3)),
            mono(((1, 4), (0, 0)), F(-1, 3)), mono(((1, 1), (0, 1))))
        source = {'modes': 8, 'particles': 4, 'hamiltonian': encode(h)}
        off, _ = coupled_fixture(source, [[0, 1], [2, 3]], F(0))
        on, receipt = coupled_fixture(source, [[0, 1], [2, 3]], F(1))
        self.assertEqual(decode(on['hamiltonian'], 8, 4), canonical(h))
        self.assertEqual(len(off['hamiltonian']), 2)
        self.assertEqual(receipt['between_terms_changing_fragment_charge'], 2)
        self.assertFalse(receipt['local_particle_number_constraints_added'])
        for partition in [[[0, 1], [1, 2, 3]], [[0, 1], [2]]]:
            with self.assertRaises(ValueError): coupled_fixture(source, partition, F(1))
        with self.assertRaises(ValueError): coupled_fixture(source, [[0, 1], [2, 3]], F(2))

    def test_endpoint_identity_and_explicit_state_continuation(self):
        from research.interacting_scaling_20260915 import budget
        from research.interacting_scaling_20260915.coupling import initialize
        from research.interacting_scaling_20260915.state import seed
        from research.molecular_collective_20260913.core import digest
        data = {'modes': 8, 'particles': 4, 'hamiltonian': encode(mono(((1, 0), (0, 0))))}
        with tempfile.TemporaryDirectory() as temporary, patch.object(budget, 'OUT', Path(temporary)/'output'):
            source = Path(temporary)/'source'; source.mkdir()
            (source/'fixture.json').write_text(json.dumps(data))
            off = initialize(source, 'off', F(0), widths=[2])
            (off/'mps').mkdir()
            off_data = json.loads((off/'fixture.json').read_text())
            old_state = seed(off_data)
            (off/'mps/state.json').write_text(json.dumps(old_state))
            on = initialize(source, 'on', F(1), widths=[2], previous=off)
            self.assertEqual(json.loads((on/'fixture.json').read_text()), data)
            new_state = json.loads((on/'continuation_seed.json').read_text())
            self.assertEqual(new_state['fixture_sha256'], digest(data))
            self.assertEqual(new_state['tensors'], old_state['tensors'])
            with self.assertRaises(ValueError): initialize(source, 'wrong_partition', F(1), 1, [2], previous=off)


if __name__ == '__main__':
    unittest.main()
