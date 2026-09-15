from fractions import Fraction as F
from itertools import combinations
import json
from math import comb
from pathlib import Path
import unittest

from experiments.marginal_dual_exact import exact_psd, functional, full_cubic_blocks, rows_for, verify
from experiments.marginal_symbolic import adj, mono, product
from tests.test_marginal_hunt_car import act


class ExactMomentTests(unittest.TestCase):
    def test_psd_handles_singular_and_indefinite_exactly(self):
        self.assertEqual(exact_psd([[1,1],[1,1]]),(True,1))
        self.assertFalse(exact_psd([[0,1],[1,1]])[0])
        self.assertFalse(exact_psd([[1,2],[2,1]])[0])
        self.assertEqual(exact_psd([[0,0],[0,1]]),(True,1))
        with self.assertRaises(ValueError):exact_psd([[1,0],[1,1]])

    def test_physical_moments_match_independent_signed_trace(self):
        modes=4;n=2;rows=rows_for(modes)
        values=[]
        for w in rows:
            k=len(w)//2;left=tuple(i for c,i in w if c);right=tuple(i for c,i in w if not c)
            values.append(F((-1)**(k*(k-1)//2)*comb(n,k),comb(modes,k)) if left==right else F(0))
        evaluate=functional(rows,values,modes)
        states=[sum(1<<i for i in indices) for indices in combinations(range(modes),n)]
        for block in full_cubic_blocks(modes):
            for u in block['words']:
                for v in block['words']:
                    p=product(adj(mono(u)),mono(v))
                    trace=sum(act(p,s).get(s,0) for s in states)/F(len(states))
                    self.assertEqual(evaluate(p),trace)

    def test_saved_witness_proves_strict_relaxation_gap(self):
        root=Path(__file__).resolve().parents[1]
        c=json.loads((root/'results/marginal_dual/witness.json').read_text())
        receipt=verify(c)
        ref=json.loads((root/'results/marginal_sector_reference/reference.json').read_text())
        self.assertGreater(F(ref['lower'])-F(receipt['energy']),F(23,1000))
        self.assertLess(F(receipt['energy']),F(328,100))

    def test_normalization_and_hamiltonian_tampering_rejected(self):
        root=Path(__file__).resolve().parents[1]
        c=json.loads((root/'results/marginal_dual/witness.json').read_text())
        c['moments'][0]='2'
        with self.assertRaises(ValueError):verify(c)
        c['moments'][0]='1';c['t']='1'
        with self.assertRaises(ValueError):verify(c)


if __name__=='__main__':unittest.main()
