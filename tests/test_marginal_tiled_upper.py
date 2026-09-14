import copy
import json
import unittest
from pathlib import Path
from fractions import Fraction as F

from experiments.marginal_tiled_upper import replay_upper8, replay_tiling, remainder_upper

UPPER8 = F(-120558410244999997006043477431315508062486123707,
           28461740687412294825487447957163629238443780483)


ROOT = Path(__file__).resolve().parents[1]


class TiledUpperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads((ROOT / 'results/marginal_graded_hubbard8/singlet_moment_sharp/certificate.json').read_text())
        cls.hamiltonian = json.loads((ROOT / 'results/marginal_graded_hubbard8/hamiltonian.json').read_text())

    def test_known_upper_replays_exactly(self):
        result = replay_upper8(self.certificate, self.hamiltonian)
        self.assertEqual(F(result['upper']), UPPER8)
        self.assertTrue(result['accepted'])

    def test_rejects_altered_hamiltonian_before_oracle(self):
        h = copy.deepcopy(self.hamiltonian)
        h['hamiltonian'][0]['coefficient'] = '5'
        with self.assertRaises(ValueError): replay_upper8(self.certificate, h)
        for field,value in [('modes',18),('particles',6)]:
            with self.assertRaises(ValueError):
                replay_upper8(self.certificate,dict(self.hamiltonian,**{field:value}))

    def test_rejects_bad_coefficients(self):
        c = copy.deepcopy(self.certificate)
        c['upper_chebyshev_coefficients'][0][0] = 10**16
        with self.assertRaises(ValueError): replay_upper8(c, self.hamiltonian)

    def test_rejects_bad_boundary_data(self):
        for mutate in (
            lambda c: c['embedding']['basis'][0].update({'0': 1}),
            lambda c: c['embedding']['basis'][0].update({'39321': 1.5}),
            lambda c: c['embedding']['basis'][0].update({'39321': 10**16}),
        ):
            c = copy.deepcopy(self.certificate); mutate(c)
            with self.assertRaises(ValueError): replay_upper8(c, self.hamiltonian)

    def test_rejects_shape_and_colliding_state_labels(self):
        c = copy.deepcopy(self.certificate)
        c['upper_chebyshev_coefficients'] = [[1]]
        with self.assertRaises(ValueError): replay_upper8(c, self.hamiltonian)
        c = copy.deepcopy(self.certificate)
        state = next(iter(c['embedding']['basis'][0]))
        c['embedding']['basis'][0][int(state)] = 1
        with self.assertRaises(ValueError): replay_upper8(c, self.hamiltonian)

    def test_tiling_remainder_requires_physical_local_certificate(self):
        with self.assertRaises(ValueError): replay_tiling(self.certificate, self.hamiltonian, 10)
        bad = {'sites': 2, 'U': '3', 't': '1', 'upper_vector': {'3': 1}}
        with self.assertRaises(ValueError): replay_tiling(self.certificate, self.hamiltonian, 10, bad)

    def test_exact_nonzero_remainder_and_extensive_upper(self):
        # In the symmetric ionic/singlet basis the two-site Hamiltonian is
        # [[4,-2],[-2,0]]. The trial coefficients (1,2) give -4/5.
        local = {'sites':2,'U':'4','t':'1','upper_vector':{3:1,6:-2,9:2,12:1}}
        self.assertEqual(F(remainder_upper(local,2)['upper']),F(-4,5))
        result = replay_tiling(self.certificate,self.hamiltonian,10,local)
        self.assertEqual(F(result['upper']),UPPER8-F(4,5))
        self.assertEqual(F(result['upper_per_site']),(UPPER8-F(4,5))/10)
        for field,value in [('upper_vector',{0:1}),('upper_vector',{3:0}),
                            ('onsite_profile',['4','4']),('upper_vector',{3:1,'3':1})]:
            with self.subTest(field=field,value=value):
                with self.assertRaises(ValueError):
                    remainder_upper(dict(local,**{field:value}),2)


if __name__ == '__main__': unittest.main()
