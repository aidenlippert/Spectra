from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_symbolic import decode, hermitian
from experiments.marginal_transfer_verify import replay


class MoleculeAcceptanceTests(unittest.TestCase):
    def test_compact_molecular_interval_and_fixture_binding(self):
        out=Path(__file__).resolve().parents[1]/'results/marginal_molecule'
        c=json.loads((out/'h4_compact_certificate.json').read_text())
        f=json.loads((out/'h4_compact_fixture.json').read_text())
        h=decode(c['hamiltonian'],8,4)
        self.assertEqual(h,decode(f['hamiltonian'],8,4))
        self.assertTrue(hermitian(h))
        self.assertTrue(all((10**12*v).denominator==1 for v in h.values()))
        result=replay(c)
        self.assertLess(F(result['width']),F(1,100000))
        reference=f['electronic_fci_total']
        self.assertLess(float(F(result['lower'])),reference)
        self.assertLessEqual(reference,float(F(result['upper']))+1e-12)

    def test_quantization_budget_covers_actual_coefficient_change(self):
        out=Path(__file__).resolve().parents[1]/'results/marginal_molecule'
        old=json.loads((out/'h4_rectangle_sto3g.json').read_text())
        new=json.loads((out/'h4_compact_fixture.json').read_text())
        a=decode(old['hamiltonian'],8,4);b=decode(new['hamiltonian'],8,4)
        error=sum(abs(a.get(w,0)-b.get(w,0)) for w in set(a)|set(b))
        self.assertEqual(error,F(new['additional_quantization_l1']))
        self.assertGreaterEqual(F(new['numerical_integral_perturbation_l1_upper']),error)


if __name__=='__main__':unittest.main()
