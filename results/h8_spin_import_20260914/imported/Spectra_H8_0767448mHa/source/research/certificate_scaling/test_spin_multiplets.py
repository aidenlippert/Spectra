import unittest
from fractions import Fraction as F
from experiments.marginal_symbolic import mono,add,scale,product,canonical,expand_squares
from experiments.marginal_hunt_car import adj
from experiments.marginal_spin_reduction import spin_operators,commutator
from research.certificate_scaling.spin_multiplets import multiplet,invariant_square,export_invariant_square
from research.certificate_scaling.test_dressed_quadratic_certificate import occupation_matrix

class SpinMultipletTests(unittest.TestCase):
    def test_all_ranks_invariant_and_exact_integer_export(self):
        for rank,weights in [(0,[1]),(1,[1,1]),(2,[2,1,2]),(3,[3,1,1,3])]:
            p=mono(tuple((1,2*i) for i in range(rank)),F(2,7))
            descendants,actual=multiplet(p,6,rank)
            self.assertEqual(actual,weights)
            sos=invariant_square(p,6,rank)
            for generator in spin_operators(6):self.assertFalse(commutator(generator,sos))
            blocks,d=export_invariant_square(p,6,rank)
            expanded,_=expand_squares(blocks,d,6)
            self.assertEqual(expanded,sos)

    def test_half_spin_creation_square_includes_hole_constant(self):
        actual=invariant_square(mono(((1,0),)),2,1)
        self.assertEqual(actual,add(mono((),2),mono(((1,0),(0,0)),-1),mono(((1,1),(0,1)),-1)))

    def test_nonorthogonal_copies_and_independent_matrix_commutation(self):
        p=add(mono(((1,0),(1,2)),F(2,3)),mono(((1,0),(1,4)),F(-3,7)))
        sos=invariant_square(p,6,2)
        _,A=occupation_matrix(sos,6,2)
        for op in spin_operators(6):
            _,B=occupation_matrix(op,6,2);n=len(A)
            self.assertTrue(all(sum(A[i][k]*B[k][j]-B[i][k]*A[k][j] for k in range(n))==0 for i in range(n) for j in range(n)))

    def test_wrong_weights_break_invariance_and_bad_highest_refuses(self):
        p=mono(((1,0),(1,2),(1,4)))
        descendants,_=multiplet(p,6,3)
        wrong=add(*(product(canonical(adj(q)),q) for q in descendants))
        self.assertTrue(commutator(spin_operators(6)[0],wrong))
        with self.assertRaises(ValueError):multiplet(mono(((1,1),)),2,1)
        with self.assertRaises(ValueError):multiplet(p,6,2)
        with self.assertRaises(ValueError):multiplet({((1,0),):1.0},2,1)

if __name__=='__main__':unittest.main()
