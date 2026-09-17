import copy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest
from research.positive_cone_20260916.cone import (one_spin_model,apply_map,integer_psd,verify,original_to_matrix_label)


class ConeTests(unittest.TestCase):
    def test_hopping_dimension_and_bipartite(self):
        labels,K,edges,eta=one_spin_model()
        self.assertEqual(len(labels),70);self.assertEqual(len(edges),10)
        self.assertTrue(all(eta[i]*eta[j]==-1 for i,j in edges))
        self.assertTrue(all(a.bit_count()==4 for a in labels))
        self.assertEqual(K,[list(row) for row in zip(*K)])

    def test_actual_signed_particle_hole_mapping_on_all_transitions(self):
        # Validation ONLY: independently instantiate physical configurations.
        from research.constructive_response_20260916.interacting_probe import rational_rows
        for rungs in (2,4):
            labels,K,edges,eta=one_spin_model(rungs)
            ix={a:i for i,a in enumerate(labels)}
            original,rows=rational_rows(rungs,F(1))
            mapped=[original_to_matrix_label(x,rungs) for x in original]
            self.assertEqual(len({(a,b) for a,b,s in mapped}),len(original))
            original_index={(a,b):i for i,(a,b,s) in enumerate(mapped)}
            for col,(a,b,phase) in enumerate(mapped):
                i,j=ix[a],ix[b];expected={}
                diagonal=8*(rungs-(a&b).bit_count())
                if diagonal:expected[col]=diagonal
                for k,c in enumerate(labels):
                    if K[k][i]:
                        row=original_index[c,b]
                        expected[row]=K[k][i]*phase*mapped[row][2]
                    if K[k][j]:
                        row=original_index[a,c]
                        expected[row]=K[k][j]*phase*mapped[row][2]
                self.assertEqual(expected,rows[col])

    def test_psd_gate_and_need_for_strict_amplitude(self):
        integer_psd([[2,1],[1,2]],strict=True)
        integer_psd([[1,1],[1,1]])
        for bad,strict in (([[1,1],[1,1]],True),([[0,1],[1,0]],False),([[1,2],[0,1]],False)):
            with self.assertRaises(ValueError):integer_psd(bad,strict)
        # Without strictness, a singular C supported only at energy +2 could
        # claim lower +2 although L(C)=KC+CK has a -2 eigenvalue.
        C=[[0,0],[0,1]];K=[[-1,0],[0,1]]
        residual=[[sum(K[i][k]*C[k][j]+C[i][k]*K[k][j] for k in range(2))-2*C[i][j]
                   for j in range(2)] for i in range(2)]
        integer_psd(residual)
        with self.assertRaises(ValueError):integer_psd(C,strict=True)

    def test_map_self_adjoint_and_adjoint_symmetry(self):
        labels,K,_,_=one_spin_model(2);n=len(labels)
        X=[[2*i-j for j in range(n)] for i in range(n)]
        Y=[[i+j*j for j in range(n)] for i in range(n)]
        LX=apply_map(X,labels,K);LY=apply_map(Y,labels,K)
        self.assertEqual(sum(X[i][j]*LY[i][j] for i in range(n) for j in range(n)),
                         sum(LX[i][j]*Y[i][j] for i in range(n) for j in range(n)))
        LT=apply_map([list(r) for r in zip(*X)],labels,K)
        self.assertEqual(LT,[list(r) for r in zip(*LX)])

    def test_accept_and_model_refusals(self):
        root=Path(__file__).resolve().parents[2]
        payload=json.loads((root/'results/positive_cone_20260916/cold/certificate.json').read_text())
        result=verify(payload)
        self.assertTrue(result['meets_target'])
        self.assertLess(F(result['width_over_t']),F(1,1000000))
        cases=[]
        p=copy.deepcopy(payload);p['U']='9';cases.append(p)
        p=copy.deepcopy(payload);p['target_spin_populations']=[3,5];cases.append(p)
        p=copy.deepcopy(payload);p['lower_over_t']='-3';cases.append(p)
        p=copy.deepcopy(payload);p['integer_C'][0][1]+=1;cases.append(p)
        p=copy.deepcopy(payload);p['integer_C']=[[0]*70 for _ in range(70)];cases.append(p)
        for p in cases:
            with self.assertRaises(ValueError):verify(p)

if __name__=='__main__':unittest.main()
