import json
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_hunt_car import mono, add
from experiments.marginal_orbital_rotation import rotation, rotate, car_check, run
from experiments.marginal_spin_reduction import spin_operators, commutator
from experiments.marginal_symbolic import decode

class OrbitalRotationTests(unittest.TestCase):
    def test_rational_rotation_inverse_recovers_one_and_two_body_polynomial(self):
        u=rotation(8); v=rotation(8,c=F(3,5),s=-F(4,5))
        h=add(mono(((1,0),(0,2)),F(7,11)), mono(((1,0),(1,4),(0,2),(0,6)),F(5,13)))
        self.assertEqual(rotate(rotate(h,u),v),h)

    def test_car_and_spin_preserving_map(self):
        self.assertTrue(car_check(rotation(12)))
        d=json.loads(Path('results/marginal_h6/direct_spin/symmetric_hamiltonian.json').read_text())
        h=decode(d['hamiltonian'],12,4); hr=rotate(h,rotation(12))
        plus,minus,z=spin_operators(12)
        self.assertFalse(commutator(hr,plus)); self.assertFalse(commutator(hr,z))

    def test_refuse_invalid_rotation(self):
        with self.assertRaises(ValueError): rotation(7)
        with self.assertRaises(ValueError): rotation(8,c=F(1,2),s=F(1,2))

    def test_actual_h6_diagnostic(self):
        r=run()
        self.assertEqual(r['sz0_dimension'],400)
        self.assertTrue(r['car_valid'])
        self.assertLess(r['spectrum_max_abs_error'],1e-10)
        self.assertIn('comparison_floor_delta',r)

if __name__=='__main__': unittest.main()
