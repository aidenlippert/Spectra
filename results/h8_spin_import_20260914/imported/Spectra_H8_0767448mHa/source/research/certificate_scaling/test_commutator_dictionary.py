"""Independent algebra checks for the H-only polynomial Gram dictionary."""
from fractions import Fraction as F
import unittest
import contextlib,io,json,tempfile
from pathlib import Path
import numpy as np
from experiments.marginal_symbolic import mono,add,scale,product,canonical,expand_squares
from experiments.marginal_hunt_car import adj
from research.certificate_scaling.commutator_dictionary import commutator,gram_columns,export_groups,polynomial_groups,diagonal_commutator,run

class CommutatorDictionaryControls(unittest.TestCase):
    def test_scs_and_construction_receipt_on_known_positive_hamiltonian(self):
        h=product(mono(((1,0),(0,0))),mono(((1,1),(0,1))))
        with tempfile.TemporaryDirectory() as root,contextlib.redirect_stdout(io.StringIO()):
            out=Path(root)/'solve';r=run(h,4,2,out,solver_seconds=2,enrich=False,solver='SCS',row_condition=True,map_backend='contraction',basis_condition='whiten')
            self.assertGreater(F(r['exact']['lower']),F(-1,10000))
            self.assertLessEqual(F(r['exact']['lower']),0)
            pre=json.loads((out/'pre_solve.json').read_text())
            self.assertEqual(pre['gram_map_nonzeros'],r['gram_map_nonzeros'])
            self.assertEqual(pre['solver'],'SCS')

    def test_export_congruence_orientation_in_original_basis(self):
        p=mono(((0,0),),F(1,3));q=product(mono(((1,1),(0,1))),mono(((0,0),)))
        W=np.array([[2.,1.],[0.,1.]])
        Q=np.array([[4.,2.],[2.,1.]]) # row [2,1], effective original row [4,3]
        blocks,den,_=export_groups([{'name':'minus','polynomials':[p,q]}],[Q],transforms=[W])
        expanded,_=expand_squares(blocks,den,2)
        factor=add(scale(p,4),scale(q,3))
        self.assertEqual(expanded,product(canonical(adj(factor)),factor))

    def test_diagonal_driver_can_add_a_new_polynomial_direction(self):
        p=add(mono(((1,0),(0,1),(0,2))),mono(((1,0),(0,1),(0,3))))
        image=diagonal_commutator(p,{0:F(0),1:F(1),2:F(2),3:F(5)})
        self.assertEqual(set(image),set(p))
        self.assertEqual(len({image[w]/p[w] for w in p}),2)

    def test_diagonal_driver_shortcut_matches_exact_car(self):
        energies={i:F(2**i,3) for i in range(4)}
        h=add(*(mono(((1,i),(0,i)),v) for i,v in energies.items()))
        p=add(mono(((1,1),(0,2),(0,3)),F(2,5)),mono(((1,0),(0,0),(0,1)),F(-3,7)))
        self.assertEqual(diagonal_commutator(p,energies),add(product(h,p),scale(product(p,h),-1)))

    def test_one_body_krylov_adds_exact_new_cubic_direction(self):
        n0=mono(((1,0),(0,0)));n1=mono(((1,1),(0,1)))
        h1=add(mono(((1,0),(0,2))),mono(((1,2),(0,0))))
        h=add(h1,product(n0,n1))
        p={w:c for w,c in commutator(h,((0,0),)).items() if len(w)==3}
        expected=add(product(h1,p),scale(product(p,h1),-1))
        def normalized(q):
            first=min(q);return {w:c/q[first] for w,c in q.items()}
        old,_=polynomial_groups(h,4,creator_channels=True)
        new,_=polynomial_groups(h,4,creator_channels=True,one_body_steps=1)
        self.assertFalse(any(normalized(q)==normalized(expected) for g in old for q in g['polynomials']))
        self.assertTrue(any(normalized(q)==normalized(expected) for g in new for q in g['polynomials']))
        self.assertTrue(all(len(w)<=3 for g in new for q in g['polynomials'] for w in q))
        with self.assertRaises(ValueError):polynomial_groups(h,4,one_body_steps=-1)

    def test_density_interaction_commutator(self):
        n0=mono(((1,0),(0,0)));n1=mono(((1,1),(0,1)))
        h=product(n0,n1);actual=commutator(h,((0,0),))
        self.assertEqual(actual,scale(product(n1,mono(((0,0),))),-1))
        self.assertTrue(all(len(w)<=3 for w in actual))

    def test_all_square_coefficients_and_export(self):
        p=mono(((0,0),),F(1,3));q=product(mono(((1,1),(0,1))),mono(((0,0),)))
        polys=[p,q];columns,indices,_=gram_columns(polys)
        reconstructed=add(*columns) # Gram = ones(2,2)
        factor=add(p,q);expected=product(canonical(adj(factor)),factor)
        self.assertEqual(reconstructed,expected)
        self.assertTrue(any(len(w)==4 for w in reconstructed))
        blocks,den,_=export_groups([{'name':'minus','polynomials':polys}],[np.ones((2,2))])
        expanded,_=expand_squares(blocks,den,2)
        self.assertEqual(expanded,expected)

    def test_cubic_square_degree_six_is_not_discarded(self):
        p=mono(((1,0),(0,1),(0,2)))
        columns,_,_=gram_columns([p])
        self.assertTrue(any(len(w)==6 for w in columns[0]))
        self.assertEqual(columns[0],product(canonical(adj(p)),p))

    def test_creator_channels_preserve_commutator_span_and_charge(self):
        ns=[mono(((1,i),(0,i))) for i in range(4)]
        h=add(product(ns[0],ns[1]),scale(product(ns[0],ns[2]),F(2,3)))
        groups,_=polynomial_groups(h,4,creator_channels=True)
        candidates=[p for g in groups for p in g['polynomials'] if all(len(w)==3 for w in p)]
        for g in groups:
            self.assertEqual(len({sum(2*c-1 for c,_ in w) for p in g['polynomials'] for w in p}),1)
        cubic={w:c for w,c in commutator(h,((0,0),)).items() if len(w)==3}
        recovered={}
        for k in {w[0][1] for w in cubic}:
            piece={w:c for w,c in cubic.items() if w[0][1]==k}
            match=next(p for p in candidates if set(p)==set(piece))
            word=next(iter(piece));recovered=add(recovered,scale(match,piece[word]/match[word]))
        self.assertEqual(recovered,cubic)

if __name__=='__main__':unittest.main()
