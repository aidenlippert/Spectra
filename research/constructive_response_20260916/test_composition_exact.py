import unittest
from fractions import Fraction as F
from research.constructive_response_20260916.composition_exact import (
    additive_error_counterexample,exact_composition,nested_exact,is_psd,
    response_bounds,schur,
)


class CompositionTests(unittest.TestCase):
    def test_noncommuting_metric_and_symmetric_response(self):
        out=exact_composition()
        self.assertNotEqual(out['response_lower'],out['response_upper'])

    def test_two_actual_eliminations_match_direct(self):
        out=nested_exact()
        self.assertEqual(out['nested'],out['direct'])
        self.assertNotEqual(out['lower'],out['upper'])

    def test_error_can_be_amplified(self):
        self.assertEqual(additive_error_counterexample(),18)

    def test_refuses_missing_positivity_premises(self):
        self.assertFalse(is_psd([[F(1),F(1)],[F(0),F(1)]]))
        self.assertFalse(is_psd([[F(0),F(1)],[F(1),F(0)]]))
        with self.assertRaises(ValueError):
            response_bounds([[F(1)]],[[F(1)]],[[F(0)]],[[F(3)]],F(1,2))
        with self.assertRaises(ValueError):
            schur([[F(1),F(0)],[F(0),F(-1)]],[0])


if __name__=='__main__':
    unittest.main()
