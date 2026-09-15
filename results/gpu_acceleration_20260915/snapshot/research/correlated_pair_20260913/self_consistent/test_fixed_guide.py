from fractions import Fraction as F
import copy
import unittest
from experiments.marginal_symbolic import mono,add,adj,canonical,product
from research.correlated_pair_20260913.self_consistent.fixed_guide import tau_map,cross,expand,ph,build,check

class FixedGuide(unittest.TestCase):
    def test_cross_unequal_weights(self):
        weights=list(map(F,[1,2,3,4]));terms=[((1,0),(1,1),(0,3),(0,2)),((1,0),(1,1),(1,2),(0,3)),((1,0),(1,1),(1,2),(1,3)),((1,0),(0,1)),((1,0),(1,1))]
        for word in terms:
            q=mono(word,F(3,7));V=canonical(add(q,adj(q)));taus=tau_map(V,weights)
            self.assertEqual(cross(taus,weights),V)
            sos,_=expand(taus,weights)
            literal=add(*(add({}) for _ in []))
            for i,t in enumerate(taus):
                from experiments.marginal_symbolic import scale
                a=add(mono(((0,i),)),t)
                literal=add(literal,scale(add(product(adj(a),a),product(t,adj(t))),weights[i]))
            self.assertEqual(canonical(literal),sos)
    def test_ph_involution(self):
        q=canonical(mono(((1,0),(1,2),(0,2),(0,0)),F(2,3)))
        self.assertEqual(ph(ph(q,2),2),q)
    def test_bad_inputs(self):
        with self.assertRaises(ValueError):tau_map({():F(1)},[F(1)])
        with self.assertRaises(ValueError):tau_map({},[F(-1)])

if __name__=='__main__':unittest.main()
