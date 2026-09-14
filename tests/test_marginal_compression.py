"""Charge splitting must preserve every word and the twirled square."""
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_compression import charge, compressed_dictionaries, hopping_components
from experiments.marginal_coefficient import dagger, dictionaries
from experiments.marginal_collective import hopping_polynomial, verify_interval
from experiments.marginal_symbolic import add, mono, product


class CompressionTests(unittest.TestCase):
    def test_components_follow_actual_asymmetric_hopping(self):
        self.assertEqual(hopping_components(hopping_polynomial(6,F(1,5)),6),[[0,3],[1,4],[2,5]])
        self.assertEqual(hopping_components(hopping_polynomial(6,F(1,5),True),6),[[0,1,3,4],[2,5]])

    def test_nonconserving_interaction_rejects_proposed_split(self):
        w=((1,0),(1,1),(0,3),(0,2))
        h=add(mono(w),mono(dagger(w)))
        with self.assertRaises(ValueError): hopping_components(h,4)

    def test_split_preserves_words_and_compresses_largest_block(self):
        h=hopping_polynomial(8,F(1,5))
        full=dictionaries(8,'mixed'); split=compressed_dictionaries(h,8,split=True)
        self.assertEqual(sorted(w for b in full for w in b['words']),sorted(w for b in split for w in b['words']))
        self.assertEqual(max(len(b['words']) for b in full),232)
        self.assertEqual(max(len(b['words']) for b in split),28)

    def test_charge_projection_equals_split_square(self):
        h=hopping_polynomial(4,F(1,5)); groups=hopping_components(h,4)
        member={i:g for g,group in enumerate(groups) for i in group}
        words=dictionaries(4,'mixed')[-2]['words']
        polynomial=add(*(mono(w,i+1) for i,w in enumerate(words)))
        adjoint=add(*(mono(dagger(w),i+1) for i,w in enumerate(words)))
        square=product(adjoint,polynomial)
        projected={w:c for w,c in square.items() if not any(charge(w,member,len(groups)))}
        pieces={}
        for i,w in enumerate(words):
            q=charge(w,member,len(groups))
            pieces[q]=add(pieces.get(q,{}),mono(w,i+1))
        split=add(*(product(add(*(mono(dagger(w),c) for w,c in p.items())),p) for p in pieces.values()))
        self.assertEqual(projected,split)

    def test_edge_and_density_dictionaries_are_adjoint_closed_subsets(self):
        h=hopping_polynomial(8,F(1,5))
        density=compressed_dictionaries(h,8,'density')
        edge=compressed_dictionaries(h,8,'edge')
        full=dictionaries(8,'mixed')
        self.assertLess(set(density[-2]['words']),set(edge[-2]['words']))
        self.assertLess(set(edge[-2]['words']),set(full[-2]['words']))
        for blocks in (density,edge):
            self.assertEqual([dagger(w) for w in blocks[-2]['words']],blocks[-1]['words'])

    def test_saved_split_bound_is_tight(self):
        root=Path(__file__).resolve().parents[1]
        cert=json.loads((root/'results/marginal_compression/m8_t1_5_full_split.json').read_text())
        self.assertLess(F(verify_interval(cert)['width']),F(1,100000))


if __name__=='__main__': unittest.main()
