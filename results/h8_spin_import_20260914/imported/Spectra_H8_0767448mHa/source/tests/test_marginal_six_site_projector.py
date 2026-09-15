import unittest
import json
from pathlib import Path
from fractions import Fraction as F
from math import comb

from experiments.marginal_projector_extendibility import overlap_grams, projector_bound, replay


class SixSiteProjectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=Path(__file__).resolve().parents[1] / 'results/marginal_graded_hubbard8'
        cls.source=json.loads((cls.base/'density_transfer/window_ceiling_certificate.json').read_text())['physical_state']
    def setUp(self):
        # A fixed-spin, half-filled six-site determinant is sufficient to
        # construct every physical embedding column independently.
        self.vector = {63: 1}

    def test_window_gram_dimensions_and_support(self):
        norm, grams, fixed = overlap_grams(self.vector, 2, 6)
        self.assertEqual(norm, 1)
        self.assertTrue(fixed)
        self.assertEqual(sum(map(len, grams.values())), 8)
        self.assertEqual(max(map(len, grams.values())), 2)
        result = projector_bound(self.vector, 2, 2, 6)
        self.assertTrue(result['accepted'])
        self.assertEqual(result['gram_dimension'], 8)
        self.assertEqual(result['support_sites'], 7)

    def test_ceiling_one_rejects_cyclic_source(self):
        with self.assertRaises(ValueError):
            projector_bound(self.source, 2, 1, 6)

    def test_actual_cyclic_two_window_gram_by_independent_tensor_embedding(self):
        v={int(s):a for s,a in self.source.items()}
        columns=[{s+4096*external:a for s,a in v.items()} for external in range(4)]
        columns += [{external+4*s:a for s,a in v.items()} for external in range(4)]
        groups={}
        for column in columns:
            counts={(sum((s>>(2*i))&1 for i in range(7)),sum((s>>(2*i+1))&1 for i in range(7))) for s in column}
            self.assertEqual(len(counts),1)
            groups.setdefault(next(iter(counts)),[]).append(column)
        expected={key:[[sum(a*right.get(s,0) for s,a in left.items()) for right in group] for left in group]
                  for key,group in groups.items()}
        norm,grams,_=overlap_grams(v,2,6)
        self.assertEqual(norm,sum(a*a for a in v.values()))
        self.assertEqual(grams,expected)
        result=projector_bound(v,4,'1084870113/500000000',6)
        self.assertEqual(result['gram_dimension'],256)
        self.assertEqual(result['maximum_psd_dimension'],36)

    def test_higher_window_dimensions(self):
        expected = {2: (8, 7), 3: (48, 8), 4: (256, 9), 5: (1280, 10)}
        for windows, (dimension, support) in expected.items():
            with self.subTest(windows=windows):
                result = projector_bound(self.vector, windows, windows, 6)
                self.assertEqual(result['gram_dimension'], dimension)
                self.assertEqual(result['support_sites'], support)
                self.assertLessEqual(result['maximum_psd_dimension'], dimension)

    def test_refusals_happen_before_gram_construction(self):
        with self.assertRaises(ValueError): overlap_grams({1: 1}, 2, 6)
        with self.assertRaises(ValueError): overlap_grams(self.vector, 1, 6)
        with self.assertRaises(ValueError): overlap_grams(self.vector, 6, 6)
        with self.assertRaises(ValueError): projector_bound(self.vector, 2, 0, 6)
        with self.assertRaises(ValueError): projector_bound(self.vector, 2, 3, 6)
        with self.assertRaises(ValueError): overlap_grams(self.vector, 2, 4)
        with self.assertRaises(ValueError): overlap_grams({63: F(1)}, 2, 6)

    def test_extended_energy_refuses_target_and_symmetry_mismatches(self):
        local=json.loads((self.base/'density_transfer/certificate.json').read_text())['targets'][0]['lower_certificate']
        c={'kind':'hubbard_projector_extension_v2','target':{'U':4,'t':1,'V':'1/2'},
           'local_window':local,'vector':self.source,'chain_sites':100,'windows':2,
           'projector_sum_ceiling':2,'penalty':'1/10','penalized_lower':-4}
        with self.assertRaisesRegex(ValueError,'target'):replay(dict(c,target={'U':3,'t':1,'V':'1/2'}))
        with self.assertRaisesRegex(ValueError,'reflection'):replay(dict(c,vector={63:1}))
        with self.assertRaisesRegex(ValueError,'spin-number'):replay(dict(c,vector={63:1,1365:1}))
        with self.assertRaisesRegex(ValueError,'Nonnegative'):replay(dict(c,penalty=-1))

    def test_five_window_spin_counts_include_both_environments(self):
        norm,grams,fixed=overlap_grams(self.source,5,6)
        expected={(up,down):5*comb(4,up-3)*comb(4,down-3)
                  for up in range(3,8) for down in range(3,8)}
        self.assertTrue(fixed)
        self.assertEqual({key:len(gram) for key,gram in grams.items()},expected)
        self.assertEqual(max(map(len,grams.values())),180)
        self.assertEqual(sum(map(len,grams.values())),1280)
        self.assertTrue(all(gram[i][i]==norm for gram in grams.values() for i in range(len(gram))))


if __name__ == '__main__':
    unittest.main()
