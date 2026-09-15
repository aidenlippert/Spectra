import copy
from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import mono,add,encode,verified_residual
from research.certificate_scaling.wedge_spectral_bound import propose
from research.collective_completion_20260914.spin_replay import check


class ResidualWedgeTests(unittest.TestCase):
    def certificate(self,constant=0):
        h=add(mono((),F(constant)),*(mono(((1,i),(0,i)),F(1,2)) for i in range(4)))
        data={'modes':4,'particles':2,'hamiltonian':encode(h)}
        core={**data,'operator_degree':3,'b':'0','denominator':1,'blocks':[],'number_multiplier':[]}
        r,_=verified_residual(core)
        cert={**data,'kind':'balanced_spin_sos_v1','alpha_multiplier':[],'core':core,
              'residual_wedge_witness':propose(r,4,2)}
        return data,cert
    def test_collective_residual_bound_and_negative_energy(self):
        for constant in (0,-3):
            data,cert=self.certificate(constant);lower=F(check(data,cert)['lower']);exact=F(1+constant)
            self.assertLessEqual(lower,exact);self.assertLess(exact-lower,F(1,10**8))
    def test_wrong_wedge_indices_refused(self):
        data,cert=self.certificate();bad=copy.deepcopy(cert)
        bad['residual_wedge_witness']['bodies'][0]['components'][0]['indices']=[[99]]
        with self.assertRaises(ValueError):check(data,bad)
    def test_averaged_singlet_residual_uses_the_checked_operator(self):
        from research.collective_completion_20260914.spin_screen import check_sector
        data,cert=self.certificate(-3)
        cert.update(kind='spin_sector_sos_v1',magnetization=0,singlet=True,
                    casimir_multiplier='0',spin_twirl=True)
        lower=F(check_sector(data,cert)['lower'])
        self.assertLessEqual(lower,-2);self.assertLess(F(-2)-lower,F(1,10**8))


if __name__=='__main__':unittest.main()
