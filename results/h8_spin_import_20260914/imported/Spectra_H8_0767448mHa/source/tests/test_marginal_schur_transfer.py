import copy
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import unittest
import numpy as np

from experiments.marginal_schur_transfer import collective_matrix,square_B,lower,replay,ldl_pivots,coupling_support_check,parity_sector_checks,resolved_complement_lower
from experiments.marginal_asymmetric_reference import model
from experiments.marginal_symbolic import encode
from experiments.marginal_transfer_verify import apply_word
from experiments.marginal_sector_reference import jacobi_data

ROOT=Path(__file__).resolve().parents[1]


class SchurTransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=json.loads((ROOT/'results/marginal_native_precision/polished/certificate.json').read_text())
        cls.certificate={'kind':'matched_spin_schur_v1','modes':10,'particles':5,'epsilon':'1/1000',
                         'lower':lower()['lower'],'hamiltonian':encode(model(F(1,1000))),
                         'independent_upper':source['independent_upper']}

    def test_projected_square_matches_all_spin_bitstrings(self):
        z=np.zeros((32,6));v=np.zeros((32,32))
        for s in range(32):
            z[s,s.bit_count()]=1
            for i in range(5):v[s^(1<<i),s]=-(i-2)
        metric=np.diag([comb(5,k) for k in range(6)])
        self.assertTrue(np.array_equal(z.T@v@z,np.zeros((6,6))))
        expected=metric@np.array(square_B(),dtype=float)
        self.assertTrue(np.array_equal(z.T@v@v@z,expected))
        self.assertEqual(square_B()[0][0],F(10))
        for matrix in (square_B(),collective_matrix()):
            weighted=[[comb(5,i)*matrix[i][j] for j in range(6)] for i in range(6)]
            self.assertEqual(weighted,list(map(list,zip(*weighted))))

    def test_complementary_spin_blocks_cover_full_fermionic_spectrum(self):
        states=[s for s in range(1024) if s.bit_count()==5];index={s:i for i,s in enumerate(states)}
        h=np.zeros((252,252))
        for word,c in model(F(0)).items():
            for j,state in enumerate(states):
                target=apply_word(word,state)
                if target:
                    dest,sign=target;h[index[dest],j]+=float(c)*sign
        expected=[]
        # j=5/2,3/2,1/2 multiplicities across all fixed-N pair-charge sectors.
        for d,multiplicity in ((0,1),(1,24),(2,75)):
            diagonal,squared=jacobi_data(10,F(1,5),d)
            block=np.diag(np.array(diagonal,dtype=float))
            off=np.sqrt(np.array(squared,dtype=float));block+=np.diag(off,1)+np.diag(off,-1)
            expected.extend(np.repeat(np.linalg.eigvalsh(block),multiplicity))
        self.assertEqual(len(expected),252)
        self.assertTrue(np.allclose(np.linalg.eigvalsh(h),np.sort(expected),atol=1e-11,rtol=0))

    def test_exact_interval_improves_native_sos_width(self):
        receipt=replay(self.certificate)
        self.assertTrue(receipt['hamiltonian_bound'])
        self.assertGreater(receipt['width_float'],0)
        self.assertLess(receipt['width_float'],1e-5)
        self.assertEqual(receipt['schur_dimension'],6)

    def test_parity_and_resolved_certificates_are_exact_and_tighter(self):
        self.assertTrue(coupling_support_check())
        self.assertEqual(parity_sector_checks(),(F(19,5),F(7,2)))
        self.assertGreater(resolved_complement_lower(),F(4))
        intervals=[]
        for kind,options in [('matched_spin_schur_parity_v1',{'parity':True}),
                             ('matched_spin_schur_resolved_v1',{'resolved':True})]:
            certificate=copy.deepcopy(self.certificate);certificate['kind']=kind
            certificate['lower']=lower(**options)['lower'];receipt=replay(certificate)
            self.assertEqual(receipt['schur_dimension'],3);intervals.append(receipt['width_float'])
        self.assertLess(intervals[0],3e-6)
        self.assertLess(intervals[1],3e-7)

    def test_spin_parity_complement_thresholds_against_all_bitstrings(self):
        j2=-1.25*np.eye(32);flip=np.zeros((32,32));h=np.zeros((32,32))
        for s in range(32):
            flip[s^31,s]=1;h[s,s]=3.75+(s.bit_count()-2.5)**2
            for i in range(5):
                h[s^(1<<i),s]=-.2
                for j in range(i+1,5):
                    t=s^((1<<i)|(1<<j)) if ((s>>i)^(s>>j))&1 else s
                    j2[t,s]+=1
        eig,u=np.linalg.eigh(j2)
        for value,dimension,threshold in ((3.75,8,4.32822),(.75,5,3.8)):
            spin=u[:,np.isclose(eig,value,atol=1e-10,rtol=0)]
            parity,rot=np.linalg.eigh(spin.T@flip@spin)
            basis=spin@rot[:,parity>.5]
            self.assertEqual(basis.shape[1],dimension)
            self.assertGreaterEqual(np.linalg.eigvalsh(basis.T@h@basis)[0],threshold-1e-12)

    def test_rejects_wrong_hamiltonian_and_failed_lower(self):
        bad=copy.deepcopy(self.certificate);bad['epsilon']='1/100'
        with self.assertRaises(ValueError):replay(bad)
        bad=copy.deepcopy(self.certificate);bad['lower']='33/10'
        with self.assertRaises(ValueError):replay(bad)

    def test_upper_witness_is_recomputed_and_validated(self):
        bad=copy.deepcopy(self.certificate);bad['independent_upper']['numerical_energy']=-1000
        self.assertEqual(replay(bad)['upper'],replay(self.certificate)['upper'])
        bad['independent_upper']['amplitudes']=[1]
        with self.assertRaises(ValueError):replay(bad)

    def test_ldl_rejects_indefinite_and_nonsymmetric_matrices(self):
        self.assertIsNone(ldl_pivots([[1,2],[2,1]]))
        with self.assertRaises(ValueError):ldl_pivots([[1,0],[1,1]])
        self.assertEqual(ldl_pivots([[2,1],[1,2]]),[F(2),F(3,2)])

    def test_search_rejects_invalid_tolerance(self):
        with self.assertRaises(ValueError):lower(tolerance=F(0))
        with self.assertRaises(ValueError):lower(F(1))


if __name__=='__main__':unittest.main()
