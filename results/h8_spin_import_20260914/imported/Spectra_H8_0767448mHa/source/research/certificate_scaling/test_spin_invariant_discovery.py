import contextlib,io,tempfile,unittest
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
import numpy as np
from experiments.marginal_symbolic import add,mono,product,expand_squares,scale,canonical
from research.certificate_scaling.spin_basis import decompose_words
from research.certificate_scaling.spin_multiplets import invariant_square
from research.certificate_scaling.spin_invariant_discovery import columns,export,run

class SpinDiscoveryTests(unittest.TestCase):
    def test_invalid_sector_and_budget_refuse_before_creating_output(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'out'
            for m,n,budget in [(3,1,1),(4,5,1),(4,-1,1),(True,1,1),(4,2,float('nan')),(4,2,float('inf')),(4,2,0)]:
                with self.assertRaises(ValueError):run({},m,n,out,seconds=budget)
                self.assertFalse(out.exists())

    def test_offdiagonal_columns_and_export_match_exact_invariant_square(self):
        words=[tuple((1,i) for i in inds) for inds in combinations(range(8),3)]
        groups,_=decompose_words(words,8);g=next(g for g in groups if g['two_spin']==3)
        k=len(g['copies']);self.assertGreater(k,1)
        cs,ix,_=columns(g);total=add(*cs)
        highest=add(*(copy[0] for copy in g['copies']))
        expected=invariant_square(highest,8,3)
        self.assertEqual(total,expected)
        blocks,d,_=export([g],[np.ones((k,k))]);actual,_=expand_squares(blocks,d,8)
        self.assertEqual(actual,expected)
        self.assertEqual(len(ix),k*(k+1)//2)

    def test_export_rational_highest_and_nonfinite_refusal(self):
        g={'copies':[[mono(((0,0),),F(1,3)),mono(((0,1),),F(-1,3))]],'weights':[1,1]}
        blocks,d,_=export([g],[np.ones((1,1))]);actual,_=expand_squares(blocks,d,2)
        expected=scale(add(mono(((1,0),(0,0))),mono(((1,1),(0,1)))),F(1,9))
        self.assertEqual(actual,expected)
        with self.assertRaises(ValueError):export([g],[np.array([[np.nan]])])

    def test_original_h_positive_control_with_small_symmetry_breaking_term(self):
        # Positive n0*n1, plus positive alpha-only number term; exact ground 0.
        h=add(product(mono(((1,0),(0,0))),mono(((1,1),(0,1)))),mono(((1,0),(0,0)),F(1,1000000)))
        with tempfile.TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
            r=run(h,4,2,Path(d)/'out',solver='CLARABEL',seconds=5)
            lo=F(r['exact']['lower']);self.assertLessEqual(lo,0);self.assertGreater(lo,F(-1,10000))

if __name__=='__main__':unittest.main()
