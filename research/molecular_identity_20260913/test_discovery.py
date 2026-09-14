"""Focused numerical integration checks; not part of exact accepting replay."""
from fractions import Fraction as F
from pathlib import Path
import contextlib
import io
import tempfile
import unittest

from experiments.marginal_symbolic import mono,add,scale,product,canonical,expand_squares
from research.certificate_scaling.commutator_dictionary import append_groups,run,gram_columns
from research.molecular_identity_20260913.coupled import augment,normalize,price


class AdditionalBlockControls(unittest.TestCase):
    def test_missing_symmetry_channel_creates_a_valid_new_block(self):
        sig=lambda w:(sum(2*f-1 for f,i in w if i%2==0),)
        base=[{'name':'linear-','polynomials':[mono(((0,0),))]}]
        # Alpha creation and two beta annihilations cannot match any linear
        # annihilator, even though the total charge remains minus one.
        p=mono(((1,0),(0,1),(0,3)))
        extras=augment(base,sig,[p])
        self.assertTrue(all('merge_into' not in g for g in extras))
        merged=append_groups(base,extras,4,sig)
        self.assertEqual(len(merged),3)
        self.assertEqual({sum(2*f-1 for f,_ in next(iter(g['polynomials'][0]))) for g in extras},{-1,1})

    def test_pricing_proposes_local_rational_operators_and_binds_h(self):
        word=((1,0),(1,1),(1,2),(0,0),(0,1),(0,2))
        raw={'modes':4,'hamiltonian':[],'rows':[[],word],'dual':[1.,-1.]}
        proposed,stats=price(raw,{},4,lambda w:(),set())
        self.assertTrue(proposed);self.assertEqual(stats['supports_priced'],1)
        for item in proposed:
            self.assertLess(item['numeric_score'],0)
            self.assertTrue(all(isinstance(v,F) for v in item['operator'].values()))
            self.assertTrue(all(len(w)==3 and sum(2*f-1 for f,_ in w)==-1 for w in item['operator']))
        with self.assertRaisesRegex(ValueError,'Hamiltonian differs'):price(raw,mono(()),4,lambda w:(),set())

    def test_invalid_additional_blocks_refused(self):
        sig=lambda w:()
        p=mono(((0,0),));base=[{'name':'linear','polynomials':[p]}]
        cases=[{'name':'empty','polynomials':[]},
               {'name':'mixed','polynomials':[p,mono(((1,0),))]},
               {'name':'outside','polynomials':[mono(((0,4),))]},
               {'name':'degree','polynomials':[mono(((1,0),(1,1),(0,0),(0,1)))]},
               {'name':'float','polynomials':[{((0,0),):0.5}]},
               {'name':'target','merge_into':'missing','polynomials':[p]},
               {'name':'wrong-charge','merge_into':'linear','polynomials':[mono(((1,0),))]}]
        for extra in cases:
            with self.subTest(name=extra['name']),self.assertRaises(ValueError):append_groups(base,[extra],4,sig)

    def test_merging_adds_cross_terms_and_preserves_input(self):
        sig=lambda w:();p=mono(((0,0),));q=mono(((1,1),(0,2),(0,3)))
        base=[{'name':'linear','polynomials':[p]}, {'name':'linear+','polynomials':[mono(((1,0),))]}]
        extras=augment(base,sig,[q]);merged=append_groups(base,extras,4,sig)
        self.assertEqual(len(base[0]['polynomials']),1)
        self.assertEqual(len(merged[0]['polynomials']),2)
        cols,_,_=gram_columns(merged[0]['polynomials'])
        self.assertTrue(cols[1]);self.assertEqual(normalize(q),normalize(scale(q,F(-2,3))))

    def test_added_cubic_projector_solver_roundtrip(self):
        a=mono(((1,0),(0,1),(0,2)))
        h=add(mono(((1,1),(0,1))),mono(((1,2),(0,2))))
        with tempfile.TemporaryDirectory() as root,contextlib.redirect_stdout(io.StringIO()):
            receipt=run(h,4,2,Path(root)/'solve',enrich=False,solver_seconds=3,
                        additional_groups=[{'name':'cubic','polynomials':[a]}])
        # N=2 can occupy orbitals 0 and 3, so the exact minimum is zero.
        self.assertLessEqual(F(receipt['exact']['lower']),0)
        self.assertGreater(F(receipt['exact']['lower']),F(-1,100000))
        self.assertEqual(receipt['additional_groups'],1)
        self.assertEqual(receipt['coefficient_max_degree'],6)


if __name__=='__main__':unittest.main()
