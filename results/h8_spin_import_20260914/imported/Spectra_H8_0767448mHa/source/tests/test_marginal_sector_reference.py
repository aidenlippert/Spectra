from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import unittest

import numpy as np

from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_sector_reference import above_lower_bound, bracket, for_certificate, reference
from experiments.marginal_cubic_probe import complete_cubic
from tests.test_marginal_hunt_car import act


class SectorReferenceTests(unittest.TestCase):
    def test_all_sector_reference_matches_independent_fermionic_matrix(self):
        for modes in (4,6,8,10):
            states=[sum(1<<i for i in indices) for indices in combinations(range(modes),modes//2)]
            lookup={s:i for i,s in enumerate(states)}
            for t in (F(0),F(1,5),F(-1,3)):
                h=hopping_polynomial(modes,t)
                matrix=np.zeros((len(states),len(states)))
                for j,state in enumerate(states):
                    for target,c in act(h,state).items():
                        matrix[lookup[target],j]=float(c)
                value=np.linalg.eigvalsh(matrix)[0]
                result=reference(modes,t,F(1,10**10))
                self.assertLessEqual(result['lower_float']-1e-12,value)
                self.assertLessEqual(value,result['upper_float']+1e-12)

    def test_pd_test_handles_singular_and_indefinite_cases(self):
        self.assertTrue(above_lower_bound([F(2),F(2)],[F(1)],F(0)))
        self.assertFalse(above_lower_bound([F(2),F(2)],[F(1)],F(1)))
        self.assertFalse(above_lower_bound([F(2),F(2)],[F(1)],F(4)))
        with self.assertRaises(ValueError): above_lower_bound([F(1)],[F(1)],F(0))
        with self.assertRaises(ValueError): bracket(10,F(1,5),0,F(0))
        with self.assertRaises(ValueError): reference(5,F(1,5))

    def test_reference_binds_actual_hamiltonian(self):
        root=Path(__file__).resolve().parents[1]
        cert=json.loads((root/'results/marginal_compression/m10_t1_5_full_split.json').read_text())
        result=for_certificate(cert)
        self.assertEqual(result['isolated_ground_sector_doubles'],0)
        self.assertLess(F(result['width']),F(1,10**12))
        cert['variational_upper']['t']='1'
        with self.assertRaises(ValueError): for_certificate(cert)

    def test_complete_cubic_includes_every_pure_triple_once(self):
        blocks=complete_cubic(hopping_polynomial(10,F(1,5)),10)
        words=[w for b in blocks if b['name'].startswith('all-triples-:') for w in b['words']]
        expected={tuple((0,i) for i in indices) for indices in combinations(range(10),3)}
        self.assertEqual(len(words),len(expected))
        self.assertEqual(set(words),expected)


if __name__=='__main__': unittest.main()
