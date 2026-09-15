from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import tempfile
import unittest
from experiments.marginal_symbolic import encode, mono
from research.interacting_scaling_20260915.fragment_control import combine, fragments


class FragmentControlTests(unittest.TestCase):
    def test_charge_convolution_covers_nonconvex_charge_energies(self):
        curves = [[F(0), F(3), F(1)], [F(2), F(-1), F(4)], [F(0), F(4), F(-2)]]
        for total in range(7):
            reference = min(sum(curves[i][q] for i, q in enumerate(choice))
                for choice in product(range(3), repeat=3) if sum(choice) == total)
            value, choices = combine(curves, total)
            self.assertEqual(value, reference)
            self.assertEqual(sum(choices), total)

    def test_coupling_and_missing_orbitals_are_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            data = {'modes': 8, 'particles': 4, 'hamiltonian': encode(mono(((1, 0), (0, 4))))}
            (path/'fixture.json').write_text(json.dumps(data))
            for receipt in [{'coupling': '0', 'partition': [[0, 1], [2, 3]]},
                {'coupling': '1', 'partition': [[0, 1], [2, 3]]},
                {'coupling': '0', 'partition': [[0, 1], [2]]}]:
                (path/'coupling.json').write_text(json.dumps(receipt))
                with self.assertRaises(ValueError): fragments(path)


if __name__ == '__main__': unittest.main()
