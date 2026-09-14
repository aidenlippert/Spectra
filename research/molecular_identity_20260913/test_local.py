"""Algebra, source binding, and refusal controls for the new local constraints."""
from fractions import Fraction as F
from itertools import product as bit_patterns
from pathlib import Path
import copy
import json
import unittest

from experiments.marginal_symbolic import mono,add,scale,product,canonical
from experiments.marginal_hunt_car import adj
from research.molecular_identity_20260913.local_blocks import local_groups
from research.molecular_identity_20260913.replay import occupation_projector,replay


class ExactLocalControls(unittest.TestCase):
    def test_all_occupation_patterns_on_nonconsecutive_modes(self):
        for support in ((0,1,2),(0,3,5),(5,0,3)):
            projectors=[occupation_projector(support,bits)[1] for bits in bit_patterns((0,1),repeat=3)]
            self.assertEqual(add(*projectors),mono(()))
            for p in projectors:self.assertEqual(product(p,p),p)
            for i,p in enumerate(projectors):
                for q in projectors[i+1:]:self.assertEqual(product(p,q),{})
        with self.assertRaises(ValueError):occupation_projector((0,0,1),(1,1,0))
        with self.assertRaises(ValueError):occupation_projector((0,2,1),(1,2,0))

    def test_local_blocks_include_all_occupation_squares(self):
        sig=lambda w:(sum(2*f-1 for f,i in w if i%2==0),)
        groups=local_groups([(0,3,5)],sig)
        squares=[product(canonical(adj(p)),p) for g in groups for p in g['polynomials']]
        for bits in bit_patterns((0,1),repeat=3):
            self.assertIn(occupation_projector((0,3,5),bits)[1],squares)
        for g in groups:
            self.assertEqual(len({(sum(2*f-1 for f,_ in w),sig(w)) for p in g['polynomials'] for w in p}),1)
        self.assertEqual(groups,local_groups([(0,3,5),(0,3,5)],sig))
        with self.assertRaises(ValueError):local_groups([(0,3,0)],sig)

    def test_offdiagonal_cross_term_is_required(self):
        p=mono(((1,0),(0,1),(0,3)));q=mono(((1,2),(0,1),(0,5)))
        whole=product(canonical(adj(add(p,q))),add(p,q))
        diagonal=add(product(canonical(adj(p)),p),product(canonical(adj(q)),q))
        self.assertNotEqual(whole,diagonal)
        cross=add(product(canonical(adj(p)),q),product(canonical(adj(q)),p))
        self.assertEqual(whole,add(diagonal,cross))
        self.assertTrue(any(len({i for _,i in w})>4 for w in cross))

    def test_exact_obstruction_and_corruption(self):
        root=Path(__file__).resolve().parents[2]
        scan=json.loads((root/'results/molecular_identity_20260913/initial_scan.json').read_text())
        raw=(root/'results/certificate_scaling/commutator_dictionary/h4_exact_dual/witness.json').read_bytes()
        receipt=replay(scan,raw)
        self.assertEqual(receipt['accepted_local_negative_directions'],32)
        self.assertLess(F(receipt['negative_occupation_value']),0)
        with self.assertRaisesRegex(ValueError,'hash mismatch'):replay(scan,raw+b' ')
        bad=copy.deepcopy(scan);bad['families'][0]['best'][0]['exact_negative_square']='0'
        with self.assertRaisesRegex(ValueError,'negative direction failed'):replay(bad,raw)


if __name__=='__main__':unittest.main()
