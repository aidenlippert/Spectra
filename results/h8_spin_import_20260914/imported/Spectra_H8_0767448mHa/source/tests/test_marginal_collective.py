"""Independent exact checks of the combinatorial upper witness."""
from fractions import Fraction as F
from itertools import product as cartesian_product
import json
from pathlib import Path
import unittest

from experiments.marginal_collective import hopping_polynomial, upper_bound, verify_interval
from tests.test_marginal_hunt_car import act


class CollectiveTests(unittest.TestCase):
    def test_combinatorial_energy_matches_signed_fermion_action(self):
        for modes in (4,6,8):
            half=modes//2
            amplitudes=[(-1)**k*(k+1) for k in range(half+1)]
            state={}
            for choices in cartesian_product((0,1),repeat=half):
                word=tuple((1,i if left else i+half) for i,left in enumerate(choices))
                for bitstring,sign in act({word:1},0).items():
                    state[bitstring]=sign*amplitudes[sum(choices)]
            norm=sum(x*x for x in state.values())
            for t in (F(0),F(1,5),F(-1,3)):
                h=hopping_polynomial(modes,t)
                energy=sum(a*coefficient*state.get(dst,0) for bitstring,a in state.items()
                           for dst,coefficient in act(h,bitstring).items())
                result=upper_bound(modes,t,amplitudes)
                self.assertEqual(F(result['upper']),F(energy,norm))
                self.assertEqual(int(result['norm']),norm)

    def test_rejects_zero_or_noninteger_witness(self):
        for amplitudes in ([0,0,0],[1,.5,1],[1,2]):
            with self.assertRaises(ValueError):
                upper_bound(4,F(1,5),amplitudes)

    def test_saved_interval_and_hamiltonian_binding(self):
        root=Path(__file__).resolve().parents[1]
        c=json.loads((root/'results/marginal_coefficient/m6_t1_5_mixed.json').read_text())
        result=verify_interval(c)
        self.assertLess(F(result['width']),F(1,10000))
        c['variational_upper']['t']='1'
        with self.assertRaises(ValueError):verify_interval(c)


if __name__=='__main__':unittest.main()
