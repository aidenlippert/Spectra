from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_symbolic import decode,encode,scale
from experiments.marginal_symmetry_transfer import replay


class AsymmetricPolishTests(unittest.TestCase):
    def test_saved_polish_reuses_positive_directions_and_tightens_bound(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_asymmetric_adapt/1_1000_penalty'
        source=json.loads((root/'certificate.json').read_text())
        certificate=json.loads((root/'polished/certificate.json').read_text())
        self.assertEqual(certificate['hamiltonian'],source['hamiltonian'])
        self.assertEqual(certificate['permutations'],source['permutations'])
        for kind in ('orbit_squares','direct_squares'):
            allowed=set()
            for item in source[kind]:
                p=decode(item['polynomial'],10,4);p=scale(p,1/max(abs(c) for c in p.values()))
                p={w:F(round(c*10**14),10**14) for w,c in p.items() if round(c*10**14)}
                allowed.add(json.dumps(encode(p)))
            for item in certificate[kind]:
                self.assertIn(json.dumps(item['polynomial']),allowed)
                self.assertGreaterEqual(F(item['weight']),0)
        receipt=replay(certificate)
        previous=json.loads((root/'receipt.json').read_text())
        self.assertLess(F(receipt['width']),F(previous['width'])/2)
        self.assertLess(receipt['orbit_squares']+receipt['direct_squares'],300)


if __name__=='__main__':unittest.main()
