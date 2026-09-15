from copy import deepcopy
from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import add, encode, mono, number_shift
from research.collective_completion_20260914.spin_screen import spin_squared
from research.interacting_scaling_20260915.dual_repair import check, eliminate, nullspace, structure
from research.interacting_scaling_20260915.singlet_trace import singlet_trace


class DualRepairTests(unittest.TestCase):
    def family(self):
        spin_plus = add(mono(((1, 0), (0, 1))), mono(((1, 2), (0, 3))))
        return {'modes': 4, 'particles': 2, 'hamiltonian': [], 'spin_defect_Ha': '0',
            'spin_average': True, 'groups': [[encode(mono(())), encode(spin_plus)]],
            'ideals': [encode(number_shift(4, 2)), encode(spin_squared(4))]}

    def test_exact_kernel_and_inconsistent_affine_constraints(self):
        self.assertEqual(nullspace([[F(1), F(1)], [F(1), F(1)]]), [[F(-1), F(1)]])
        with self.assertRaises(ValueError): eliminate([({}, F(1))])
        with self.assertRaises(ValueError): eliminate([(mono(()), F(1)), (mono(()), F(2))])

    def test_positive_singlet_trace_and_acceptance_refusals(self):
        family = self.family()
        words = structure(family)[-1]
        witness = {'moments': [{'word': w, 'value': str(singlet_trace(mono(w), 4, 2))} for w in words]}
        self.assertEqual(check(family, witness)['Gram_checks'][0]['rank'], 1)
        for word, bad in [((), '2'), (((1, 0), (0, 0)), '1/7')]:
            altered = deepcopy(witness)
            for item in altered['moments']:
                if item['word'] == word: item['value'] = bad
            with self.assertRaises(ValueError): check(family, altered)
        invalid = deepcopy(family); invalid['spin_defect_Ha'] = '1/100'
        with self.assertRaises(ValueError): check(invalid, witness)


if __name__ == '__main__':
    unittest.main()
