from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from experiments.marginal_symbolic import add, decode, encode, mono, product, scale
from research.collective_completion_20260914.spin_replay import alpha_shift
from research.interacting_scaling_20260915.solve import magnetic_certificate
from research.interacting_scaling_20260915.singlet_trace import magnetic_trace
from research.interacting_scaling_20260915 import nonsinglet


class MagneticTests(unittest.TestCase):
    def certificate(self):
        h = mono(((1, 0), (0, 0)), F(3, 7))
        y = mono(((1, 2), (0, 2)), F(2, 5))
        return {'modes': 8, 'particles': 4, 'alpha_multiplier': encode(y),
            'casimir_multiplier': '0', 'spin_ladder_multiplier': [],
            'core': {'hamiltonian': encode(add(h, scale(product(alpha_shift(8, 4), y), -1)))}}, h, y

    def test_sector_shift_matches_actual_magnetic_identity(self):
        certificate, h, y = self.certificate()
        result = magnetic_certificate(certificate)
        z = add(alpha_shift(8, 4), mono((), -1))
        self.assertEqual(decode(result['core']['hamiltonian'], 8, 4), add(h, scale(product(z, y), -1)))
        self.assertEqual(magnetic_trace(product(z, y), 8, 3, 1), 0)
        self.assertFalse(result['singlet'])
        self.assertFalse(result['spin_twirl'])
        self.assertEqual(result['magnetization'], 1)

    def test_singlet_ideals_are_refused(self):
        certificate, _, _ = self.certificate()
        for key, value in [('casimir_multiplier', '1/10'),
                           ('spin_ladder_multiplier', encode(mono(((1, 1), (0, 0)))) )]:
            invalid = deepcopy(certificate)
            invalid[key] = value
            with self.assertRaises(ValueError): magnetic_certificate(invalid)

    def test_direct_magnetic_hierarchy_keeps_global_charge(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root/'source'; source.mkdir()
            (source/'fixture.json').write_text(json.dumps({'modes': 10, 'particles': 6}))
            (source/'upper.json').write_text('{}')
            with patch.object(nonsinglet, 'OUT', root):
                designs = [json.loads((nonsinglet.initialize(source, f'level{i}', i)/'design.json').read_text()) for i in range(3)]
            self.assertFalse(designs[0]['clusters'])
            self.assertEqual(designs[1]['clusters'], designs[2]['clusters'])
            self.assertFalse(designs[1]['collective_pairs'])
            self.assertTrue(designs[2]['collective_pairs'])
            self.assertTrue(all(d['global_particle_number_only'] and d['magnetization'] == 1 for d in designs))

    def test_existing_screen_cannot_masquerade_as_fresh_discovery(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary)
            (source/'nonsinglet.json').write_text('{}')
            with self.assertRaises(ValueError): nonsinglet.construct(source, 0)


if __name__ == '__main__':
    unittest.main()
