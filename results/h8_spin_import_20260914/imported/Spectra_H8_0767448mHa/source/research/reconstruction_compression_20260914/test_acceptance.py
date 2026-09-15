import copy
from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import encode,mono
from research.reconstruction_compression_20260914.replay import check

class LowerAcceptanceTests(unittest.TestCase):
    def fixture(self):
        data={'modes':6,'particles':3,'hamiltonian':encode(mono((),F(-2)))}
        cert={'modes':6,'particles':3,'hamiltonian':data['hamiltonian'],'operator_degree':3,
              'b':'-2','number_multiplier':[],'denominator':1,'blocks':[]}
        return data,cert

    def test_scalar_and_endpoint_mutation(self):
        data,cert=self.fixture();self.assertEqual(F(check(data,cert)['lower']),-2)
        cert['b']='0';rec=check(data,cert)
        self.assertEqual(F(rec['lower']),-2)
        self.assertEqual(F(rec['residual_l1']),2)

    def test_sextic_remainder_is_included(self):
        data,cert=self.fixture();cert['blocks']=[{'words':[[[0,0],[0,1],[0,2]]],'factor':[[1]]}]
        rec=check(data,cert)
        self.assertEqual(rec['residual_max_degree'],6)
        self.assertEqual(F(rec['lower']),-3)

    def test_binding_and_schema_refusals(self):
        data,base=self.fixture()
        for mutation in ('sector','H','degree','coefficient','factor','multiplier'):
            cert=copy.deepcopy(base)
            if mutation=='sector':cert['particles']=2
            if mutation=='H':cert['hamiltonian']=encode(mono((),F(-3)))
            if mutation=='degree':cert['operator_degree']=4
            if mutation=='coefficient':cert['b']=-2.
            if mutation=='factor':cert['blocks']=[{'words':[[[0,0]]],'factor':[[1.0]]}]
            if mutation=='multiplier':cert['number_multiplier']=encode(mono(((1,0),(0,1))))
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):check(data,cert)

if __name__=='__main__':unittest.main()
