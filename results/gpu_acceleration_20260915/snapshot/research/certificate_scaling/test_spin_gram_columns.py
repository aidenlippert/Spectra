import json
import unittest
from pathlib import Path
from fractions import Fraction as F
from experiments.marginal_symbolic import decode,canonical,adj,product,scale,add,mono
from experiments.marginal_coefficient import dictionaries
from research.certificate_scaling.spin_basis import decompose_words
from research.certificate_scaling.spin_parity import spin_parity
from research.certificate_scaling.spin_gram_columns import columns,adjoint


def reference(group):
    copies=group['copies'];cols=[];ix=[];work=0
    for i,left in enumerate(copies):
        for j in range(i,len(copies)):
            result={}
            for p,q,weight in zip(left,copies[j],group['weights']):
                result=add(result,scale(product(canonical(adj(p)),q),weight));work+=len(p)*len(q)
            if i!=j:result=add(result,canonical(adj(result)))
            cols.append(result);ix.append(i*len(copies)+j)
    return cols,ix,work


class SpinColumnsTests(unittest.TestCase):
    def test_every_h4_column_matches_independent_car_reference(self):
        p=Path(__file__).resolve().parents[2]/'results/certificate_scaling/active_space_ladder/h4/fixture.json'
        f=json.loads(p.read_text());m=f['modes'];masks=spin_parity(decode(f['hamiltonian'],m,4),m)
        for family in dictionaries(m,'mixed'):
            words=[w for w in family['words'] if len(w)!=1]
            if not words:continue
            groups,_=decompose_words(words,m,masks)
            for g in groups:self.assertEqual(columns(g),reference(g))

    def test_rational_crosscopy_and_adjoint_convention(self):
        p=add(mono(((1,0),(0,1),(0,2)),F(1,3)),mono(((0,3),),F(-2,7)))
        q=add(mono(((1,2),(0,0),(0,3)),F(-5,11)),mono(((0,1),),F(7,13)))
        g={'copies':[[p],[q]],'weights':[3]}
        self.assertEqual(columns(g),reference(g))
        for poly in (p,q,*reference(g)[0]):self.assertEqual(adjoint(poly),canonical(adj(poly)))
        with self.assertRaises(ValueError):adjoint(mono(((0,0),(1,0))))
        with self.assertRaises(ValueError):columns({'copies':[[p]],'weights':[1,2]})

if __name__=='__main__':unittest.main()
