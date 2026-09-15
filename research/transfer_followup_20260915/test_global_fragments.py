from fractions import Fraction as F
import unittest
from experiments.marginal_symbolic import encode
from research.transfer_followup_20260915.global_fragments import check_line, disjoint
from research.transfer_solver_20260915.independent_oracle import matrix


class GlobalChargeComposition(unittest.TestCase):
    def test_supporting_line_controls_all_composite_charge_assignments(self):
        # One two-mode fragment has exact charge energies 0, -2, -1.
        h = {((1, 0), (0, 0)): F(-2), ((1, 1), (0, 1)): F(-2),
             ((1, 0), (1, 1), (0, 1), (0, 0)): F(3)}
        fragment = {'modes': 2, 'particles': 1, 'hamiltonian': encode(h)}
        check_line({0: F(0), 1: F(-2), 2: F(-1)}, 1, F(-2), F(-1, 2), 2)
        total = {**h, **{tuple((a, i+2) for a, i in w): v for w, v in h.items()}}
        composite = {'modes': 4, 'particles': 2, 'hamiltonian': encode(total)}
        disjoint(composite, fragment, 2)
        states, H = matrix(composite)
        self.assertEqual(min(H[i][i] for i in range(len(H))), -4)
        self.assertTrue(all(not H[i][j] for i in range(len(H)) for j in range(len(H)) if i != j))
        self.assertIn(3, states)  # Both particles on one fragment were included.

    def test_charge_instability_and_missing_sector_refused(self):
        for lowers in ({0: F(-4), 1: F(-2), 2: F(-4)}, {1: F(-2), 2: F(-1)}):
            with self.assertRaises(ValueError):
                check_line(lowers, 1, F(-2), F(0), 2)


if __name__ == '__main__':
    unittest.main()
