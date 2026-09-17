from fractions import Fraction as F
import unittest
from research.constructive_response_20260916.composition_exact import (
    add,sc,mm,tr,eye,is_psd,
)
from research.constructive_response_20260916.metric_certificate import energy_bound,square_model,sharpness_counterexample


class MetricBoundTests(unittest.TestCase):
    def test_noncommuting_large_response_exact_lower(self):
        D=[[F(2),F(0)],[F(0),F(3)]]
        X=[[F(100),F(1)],[F(2),F(50)]];rho=F(1,10)
        R=sc(rho,X);B=add(mm(D,X),R)
        A=add(mm(mm(tr(X),D),X),sc(2*rho,mm(tr(X),X)))
        eta,K,M,_=energy_bound(A,B,D,X,F(2),rho)
        H=[A[i]+tr(B)[i] for i in range(2)]+[B[i]+D[i] for i in range(2)]
        self.assertTrue(is_psd(add(H,sc(eta,eye(4)))))
        self.assertEqual(eta,F(1,180))
        self.assertNotEqual(mm(D,X),mm(X,D))
        with self.assertRaises(ValueError):
            energy_bound(A,B,D,X,F(2),F(0))

    def test_refuses_invalid_gap_and_retained_sign(self):
        with self.assertRaises(ValueError):
            energy_bound([[F(0)]],[[F(0)]],[[F(1)]],[[F(0)]],F(2),F(0))
        with self.assertRaises(ValueError):
            energy_bound([[F(-1)]],[[F(0)]],[[F(2)]],[[F(0)]],F(2),F(0))

    def test_square_has_correct_particle_and_hopping_structure(self):
        labels,H,P,Q=square_model()
        self.assertEqual((len(labels),len(P),len(Q)),(36,6,30))
        self.assertTrue(all(x.bit_count()==4 for x in labels))
        self.assertTrue(all(H[i][j]==0 for i in P for j in P))

    def test_exact_optimality_and_gap_threshold(self):
        # eta*=1/180; every listed smaller proposed allowance is invalid.
        for fraction in (F(0),F(1,2),F(99,100),F(9999,10000)):
            c=sharpness_counterexample(F(2),F(1,10),fraction/F(180))
            self.assertLess(F(c['determinant_at_proposed_lower']),0)
        with self.assertRaises(ValueError):
            sharpness_counterexample(F(2),F(1,10),F(1,180))
        for delta in (F(1),F(2)):
            for attempted in (F(0),F(10),F(10**6)):
                c=sharpness_counterexample(delta,F(1),attempted)
                self.assertLess(F(c['determinant_at_proposed_lower']),0)


if __name__=='__main__':
    unittest.main()
