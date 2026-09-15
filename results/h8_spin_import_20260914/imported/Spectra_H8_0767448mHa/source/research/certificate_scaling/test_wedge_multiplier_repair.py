"""Controls for the sparse fixed-N wedge multiplier repair LP."""
import unittest
from fractions import Fraction as F

from experiments.marginal_symbolic import add, decode, encode
from research.certificate_scaling.fermionic_ratio_chain import hamiltonian, compile_h
from research.certificate_scaling.wedge_multiplier_repair import repair


class WedgeRepairControls(unittest.TestCase):
    def test_zero_ratio_control_and_hopping_mutation(self):
        h = hamiltonian(4, F(1))
        cert, _ = compile_h(h, 4, 2)
        out, rec = repair(cert, seconds=10)
        self.assertEqual(rec["after"]["lower"], "0")
        self.assertEqual(out["hamiltonian"], cert["hamiltonian"])
        mutated = dict(cert)
        mutated["number_multiplier"] = encode(add(decode(cert["number_multiplier"], 4, 4), {
            ((1, 0), (0, 1)): F(1, 10),
            ((1, 1), (0, 0)): F(1, 10),
        }))
        repaired, rec2 = repair(mutated, seconds=10)
        self.assertLess(F(rec2["before"]["lower"]), F(0))
        self.assertEqual(F(rec2["after"]["lower"]), F(0))
        self.assertEqual(repaired["hamiltonian"], cert["hamiltonian"])
        self.assertEqual(repaired["blocks"], mutated["blocks"])
        self.assertEqual(repaired["b"], mutated["b"])
        self.assertFalse(rec2["factors_changed"])


if __name__ == '__main__':
    unittest.main()
