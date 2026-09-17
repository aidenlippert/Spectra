import unittest
from fractions import Fraction as F

from research.constructive_response_20260916.dimer_family import (
    RAD, bare_weight_interval, diagnostics, dimer_energy_interval,
    dressed_parameters, dressed_weight_interval, run, sqrt_interval,
)


class DimerFamilyTests(unittest.TestCase):
    def test_local_energy_and_weight_enclosures(self):
        elo, ehi = dimer_energy_interval()
        self.assertLess(elo, F(-472, 1000))
        self.assertGreater(ehi, F(-473, 1000))
        self.assertLess(elo, ehi)
        wlo, whi = bare_weight_interval()
        self.assertLessEqual(wlo, whi)
        self.assertLess(wlo, F(1))
        self.assertLess(F(0), wlo)

    def test_normalized_rational_dressing(self):
        lo, hi = dressed_weight_interval()
        self.assertLessEqual(lo, hi)
        self.assertGreater(lo, F(0)); self.assertLess(hi, F(1))
        c, s = dressed_parameters()
        self.assertEqual(c*c + s*s, 1)

    def test_independent_eigen_equation(self):
        elo, ehi = dimer_energy_interval()
        # det(H-EI)=E^2-8E-4; the negative root is enclosed.
        self.assertGreaterEqual(elo*elo - 8*elo - 4, 0)
        self.assertLessEqual(ehi*ehi - 8*ehi - 4, 0)
        self.assertLess(elo, ehi)

    def test_dressed_expectation_is_used(self):
        bare = diagnostics(8, False)['global_Qgap_upper_over_t']
        dressed = diagnostics(8, True)['global_Qgap_upper_over_t']
        self.assertNotEqual(bare, dressed)

    def test_negative_dressed_coefficient_preserves_interval_direction(self):
        lo,hi=dressed_weight_interval(F(0),F(1))
        wlo,whi=bare_weight_interval()
        self.assertEqual((lo,hi),(1-whi,1-wlo))
        with self.assertRaises(ValueError):
            dressed_weight_interval(F(1),None)

    def test_monotone_global_overlap_and_no_enumeration(self):
        rows = [diagnostics(L) for L in (1, 2, 4, 8, 16)]
        self.assertTrue(all(rows[i]['global_overlap_interval'][1] != rows[i+1]['global_overlap_interval'][1]
                            for i in range(4)))
        self.assertTrue(all(r['global_states_constructed'] == 0 for r in rows))
        self.assertGreater(F(rows[-1]['response_norm_lower']), F(1))

    def test_sqrt_rejects_bad_inputs_and_preserves_direction(self):
        with self.assertRaises(ValueError): sqrt_interval(F(-1))
        with self.assertRaises(ValueError): sqrt_interval(F(2), 0)
        lo, hi = sqrt_interval(RAD)
        self.assertLessEqual(lo * lo, RAD)
        self.assertGreaterEqual(hi * hi, RAD)

    def test_table_shape(self):
        result = run()
        self.assertEqual([len(x['rows']) for x in result['families']], [9, 9])
        self.assertEqual(result['families'][0]['rows'][-1]['L'], 256)


if __name__ == '__main__':
    unittest.main()
