"""Small explicit oracles are tests only, never compact-upper acceptance."""
from fractions import Fraction as F
from math import lcm
import json
import unittest
from experiments.marginal_symbolic import decode
from research.global_response_20260913 import slater
from research.compact_response_20260913 import program
from research.certificate_scaling.streaming_reference_upper import upper
from research.global_response_20260913.lifted_upper.lift import lift, check as lift_check


def determinant_columns(s, na, nb, aa, bb):
    state = {0: F(1)}
    for spin, occupied, ts in ((0, na, aa), (1, nb, bb)):
        for col in range(occupied):
            i = col % (s//2); j = i+s//2; t = F(ts[i])
            coeffs = ((2*i+spin, F(1)), (2*j+spin, t)) if col < s//2 else ((2*i+spin, -t), (2*j+spin, F(1)))
            new = {}
            for bits, amp in state.items():
                for mode, coeff in coeffs:
                    if bits >> mode & 1: continue
                    target = bits | 1 << mode
                    new[target] = new.get(target, F(0)) + amp*coeff*(-1)**((bits & ((1 << mode)-1)).bit_count())
            state = {bits: a for bits, a in new.items() if a}
    return state


def oracle_upper(data, vector):
    den = lcm(*(x.denominator for x in vector.values()))
    witness = {'states': sorted(vector), 'amplitudes': [int(vector[s]*den) for s in sorted(vector)]}
    return upper(data, witness)[0]


class CompactExactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((program.ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json').read_text())

    def test_rotated_energy_and_norm_against_independent_expansion(self):
        aa, bb = [F(1, 3), F(-2, 5)], [F(-1, 4), F(3, 7)]
        v = determinant_columns(4, 2, 2, aa, bb)
        result = slater.evaluate(self.data, aa, bb, 2, 2)
        self.assertEqual(result['norm'], sum(x*x for x in v.values()))
        self.assertEqual(result['energy'], oracle_upper(self.data, v))
        self.assertEqual(len(v), 16)

    def test_both_occupied_pair_columns(self):
        r, norm = slater.projector(4, 3, [F(2, 3), F(1, 7)])
        self.assertEqual(r[0][0], 1)
        self.assertEqual(r[2][2], 1)
        self.assertEqual(r[0][2], 0)
        self.assertEqual(norm, F(13, 9)**2*F(50, 49))

    def test_nonnormal_car_and_pure_spin_refusal(self):
        r, _ = slater.projector(4, 2, [F(1, 3), F(2, 5)])
        self.assertEqual(slater.wick(((0, 0), (1, 0)), r), 1-r[0][0])
        with self.assertRaisesRegex(ValueError, 'spin'):
            slater.make_certificate(self.data, ['1/3', '0'], ['0', '0'], 2, 2, 0)
        c = slater.make_certificate(self.data, ['1/3', '-2/5'], ['1/3', '-2/5'], 2, 2, 0)
        self.assertEqual(slater.check(self.data, c)['S2'], '0')

    def test_response_lift_exact_moments_against_compiled_car(self):
        c = lift(self.data, F(7, 10))
        v = {c['hf_state']: F(1)}
        for bits, amp in zip(c['chi_support'], c['chi_coefficients']):
            v[bits] = -F(c['alpha'])*F(amp)
        self.assertGreater(c['chi_terms'], 0)
        self.assertEqual(F(c['denominator']), sum(x*x for x in v.values()))
        self.assertEqual(F(lift_check(self.data, c)['upper_Ha']), oracle_upper(self.data, v))
        self.assertEqual(c['source_term_checks'], len(decode(self.data['hamiltonian'], 8, 4))*(1+c['chi_terms']))

    def test_wrong_particle_sector_refused(self):
        bad = dict(self.data); bad['particles'] = 3
        with self.assertRaises(ValueError): lift(bad, F(1))


if __name__ == '__main__': unittest.main()
