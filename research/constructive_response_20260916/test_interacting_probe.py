import unittest
from fractions import Fraction as F
import numpy as np
from research.constructive_response_20260916.interacting_probe import (
    balanced_basis,doublons,hop,local_rotation,lifted_gate,numerical,rational_rows,
)


class HubbardDiagnosticTests(unittest.TestCase):
    def test_hopping_sign_and_adjoint(self):
        self.assertEqual(hop(0b101,1,2),(0b011,1))
        self.assertEqual(hop(0b101,3,0),(0b1100,-1))
        for label in range(16):
            for i in range(4):
                for j in range(4):
                    moved=hop(label,i,j)
                    if moved is not None:
                        back=hop(moved[0],j,i)
                        self.assertEqual(back,(label,moved[1]))

    def test_dimer_spectrum_and_charge(self):
        labels,rows=rational_rows(1,F(0))
        self.assertEqual(labels,[3,6,9,12])
        vals=np.linalg.eigvalsh(numerical(rows).toarray())
        expected=[4-2*np.sqrt(5),0,8,4+2*np.sqrt(5)]
        np.testing.assert_allclose(vals,expected,atol=1e-12)
        for label in labels:
            self.assertEqual(label.bit_count(),2)
        h=numerical(rows).toarray()
        singlet=np.array([0,-1,1,0])/np.sqrt(2)
        ionic=np.array([1,0,0,1])/np.sqrt(2)
        self.assertAlmostEqual(ionic@h@singlet,-2,places=12)

    def test_gate_is_exact_and_uses_correct_rotation_sign(self):
        g=local_rotation()
        for i in range(16):
            for j in range(16):
                self.assertEqual(sum(g[k][i]*g[k][j] for k in range(16)),int(i==j))
        labels,rows=rational_rows(1,F(0))
        G=lifted_gate(labels,0,g).toarray()
        h=numerical(rows).toarray()
        rotated=G.T@h@G
        np.testing.assert_allclose(np.linalg.eigvalsh(rotated),np.linalg.eigvalsh(h),atol=1e-12)
        spin=[i for i,x in enumerate(labels) if doublons(x,2)==0]
        charge=[i for i,x in enumerate(labels) if doublons(x,2)>0]
        self.assertLess(np.linalg.norm(rotated[np.ix_(charge,spin)]),1e-4)

    def test_global_particle_sector_includes_charge_transfer(self):
        labels,rows=rational_rows(2,F(1))
        self.assertEqual(len(labels),36)
        self.assertTrue(any((x&15).bit_count()==1 for x in labels))
        self.assertTrue(any((x&15).bit_count()==3 for x in labels))
        spin=[i for i,x in enumerate(labels) if doublons(x,4)==0]
        self.assertEqual(len(spin),6)
        for i in spin:
            for j in spin:
                self.assertEqual(rows[i].get(j,F(0)),0)
        self.assertEqual(len(balanced_basis(6)),400)

    def test_independent_exact_certificate_model_matches(self):
        from research.constructive_response_20260916.metric_certificate import square_model
        labels,rows=rational_rows(2,F(1))
        other,H,P,Q=square_model()
        self.assertEqual(labels,other)
        self.assertEqual([[row.get(j,F(0)) for j in range(len(labels))] for row in rows],H)


if __name__=='__main__':
    unittest.main()
