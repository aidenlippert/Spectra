import copy, json, tempfile, unittest
from pathlib import Path
from fractions import Fraction as F
from experiments.marginal_comparison_witness import witness_matrix, ceiling, principal_psd, replay, construct
from tests.test_marginal_localized_basis import fixture

class TestComparisonWitness(unittest.TestCase):
    def test_all_principal_psd_and_trace_formula(self):
        h=[[F(2),F(-1),F(3),F(0),F(-2)],[F(-1),F(4),F(2),F(7),F(1)],
           [F(3),F(2),F(5),F(-3),F(4)],[F(0),F(7),F(-3),F(3),F(-1)],
           [F(-2),F(1),F(4),F(-1),F(6)]]
        v=[2,3,5,7,11]
        for k in (2,3,4):
            x=witness_matrix(h,v,k)
            self.assertTrue(principal_psd(x,k))
            self.assertEqual(ceiling(h,v,k),sum(h[i][j]*x[j][i] for i in range(5) for j in range(5))/sum(z*z for z in v))

    def test_reconstructs_physical_matrix_and_exclusion_direction(self):
        data,_=fixture()
        with tempfile.TemporaryDirectory() as z:
            root=Path(z);source=root/'h.json';source.write_text(json.dumps(data))
            result=construct(source,[3],root/'out',2,F(10))
            self.assertTrue(result['target_excluded'])
            c=json.loads((root/'out/certificate.json').read_text())
            self.assertEqual(replay(c),result)
            changed=dict(c,target_lower='-100');self.assertFalse(replay(changed)['target_excluded'])
            bad=dict(c,ceiling='100');
            with self.assertRaises(ValueError): replay(bad)
            bad=copy.deepcopy(c);bad['witness_states'][0]=3
            with self.assertRaises(ValueError): replay(bad)
            bad=copy.deepcopy(c);bad['v'][0]+=1000000
            with self.assertRaises(ValueError): replay(bad)
            with self.assertRaises(ValueError): construct(source,[3],root/'out',2,F(10))

    def test_strict_vector_and_matrix_gates(self):
        for v in ([0],[True],[-1],[10**12+1]):
            with self.assertRaises(ValueError): witness_matrix([[F(1)]],v,2)
        with self.assertRaises(ValueError): witness_matrix([[F(1)]],[1],5)
        with self.assertRaises(ValueError): witness_matrix([[1,2],[0,1]],[1,1],2)
        with self.assertRaises(ValueError): replay({'kind':'comparison_family_dual_v1'})

if __name__=='__main__': unittest.main()
