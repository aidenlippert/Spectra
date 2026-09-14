from copy import deepcopy
from fractions import Fraction as F
import json
import unittest

from research.ch2_validation_20260913 import certify
from research.certificate_scaling.streaming_reference_upper import upper
from research.certificate_scaling.commutator_dual_witness import psd


class CertifiedSpinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((certify.OUT/'fixture.json').read_text())
        cls.certs=[json.loads((certify.OUT/f'spin_{s}_certificate.json').read_text()) for s in (0,1)]

    def test_two_electron_spin_decomposition_exact(self):
        S,den,_=certify.matrix(certify.spin_square(4),[3,6,9,12])
        rational=[[F(x,den) for x in row] for row in S]
        self.assertEqual(psd(rational)['rank'],1)
        self.assertEqual(sum(rational[i][i] for i in range(4)),2)
        self.assertEqual([[sum(rational[i][k]*rational[k][j] for k in range(4)) for j in range(4)] for i in range(4)],[[2*x for x in row] for row in rational])

    def test_both_model_intervals_and_pure_spin(self):
        results=[certify.check(self.data,c) for c in self.certs]
        for r in results:
            self.assertTrue(r['exact_spin_purity']);self.assertTrue(r['exact_H_commutes_S2'])
            self.assertLess(F(r['width_Ha']),F(1,10000))
        self.assertLess(F(results[1]['upper_Ha']),F(results[0]['lower_Ha']))

    def test_core_potential_and_constant_match_RHF(self):
        generation=json.loads((certify.OUT/'generation.json').read_text())
        hf,_=upper(self.data,{'states':[63],'amplitudes':[1]})
        self.assertLess(abs(hf-F(str(generation['RHF_total_Ha']))),F(self.data['operator_rounding_envelope_Ha'])+F(1,10**10))
        self.assertEqual(self.data['particles']+2*self.data['frozen_core_orbitals'],self.data['physical_electrons'])

    def test_wrong_binding_and_accounting_refused(self):
        with self.assertRaises(ValueError):certify.check(self.data,self.certs[0]|{'fixture_sha256':'wrong'})
        with self.assertRaises(ValueError):certify.prepare(self.data|{'physical_electrons':6},0)
        with self.assertRaises(ValueError):certify.check(self.data,self.certs[0]|{'factor_denominator':0})

    def test_factor_and_spin_tampering_refused(self):
        c=deepcopy(self.certs[0]);c['lower_factor'][0][0]+=10**15
        with self.assertRaises(ValueError):certify.check(self.data,c)
        c=deepcopy(self.certs[0]);c['independent_upper']['amplitudes'][0]+=1
        value,_=upper(self.data,c['independent_upper']);c['upper_Ha']=str(value)
        # Choose a coefficient in a genuinely coupled spin pattern, if the
        # first state is closed shell and its coefficient alone stays pure.
        s=c['independent_upper']['states'][0]
        if all(bool(s&(1<<(2*i)))==bool(s&(1<<(2*i+1))) for i in range(6)):
            c=deepcopy(self.certs[0])
            j=next(j for j,s in enumerate(c['independent_upper']['states']) if any(bool(s&(1<<(2*i)))!=bool(s&(1<<(2*i+1))) for i in range(6)))
            c['independent_upper']['amplitudes'][j]+=1;value,_=upper(self.data,c['independent_upper']);c['upper_Ha']=str(value)
        with self.assertRaisesRegex(ValueError,'spin pure'):certify.check(self.data,c)


if __name__=='__main__':unittest.main()
