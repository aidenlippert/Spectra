"""Independent controls for the biased fermionic ratio-chain compiler."""
from __future__ import annotations

import itertools
import unittest
from fractions import Fraction as F

from experiments.marginal_transfer_verify import apply_word
from research.certificate_scaling.fermionic_ratio_chain import (
    compile_h, hamiltonian, structured_replay,
)


def apply_h(h, state, amplitudes):
    out = {s: F(0) for s in amplitudes}
    for word, coefficient in h.items():
        for source, amplitude in amplitudes.items():
            result = apply_word(word, source)
            if result:
                target, sign = result
                out[target] += coefficient * amplitude * sign
    return out


class FermionicRatioControls(unittest.TestCase):
    def test_exact_zero_vector_for_ratios_and_all_small_sectors(self):
        for r in (F(1), F(4), F(1, 4)):
            h = hamiltonian(4, r)
            for n in range(5):
                amps = {s: r ** sum(i for i in range(4) if (s >> i) & 1)
                        for s in range(16) if s.bit_count() == n}
                image = apply_h(h, None, amps)
                self.assertTrue(all(value == 0 for value in image.values()))
                cert, rec = compile_h(h, 4, n)
                self.assertEqual(rec["lower"], "0")
                self.assertEqual(structured_replay(cert)["width"], "0")

    def test_noncommuting_adjacent_terms_and_mutation_rejection(self):
        h = hamiltonian(4, F(4))
        left = {w: v for w, v in h.items() if all(i < 2 for _, i in w)}
        right = {tuple((c, i + 1) for c, i in w): v
                 for w, v in left.items()}
        from experiments.marginal_symbolic import add, product, scale
        comm = add(product(left, right), scale(product(right, left), -1))
        self.assertTrue(comm)
        bad = dict(h); bad[((1, 0), (0, 0))] += 1
        with self.assertRaises(ValueError):
            compile_h(bad, 4, 2)

    def test_wrong_upper_ratio_and_malformed_sector_rejected(self):
        cert, _ = compile_h(hamiltonian(4, F(4)), 4, 2)
        cert["structured_upper"]["r"] = "1/4"
        with self.assertRaises(ValueError):
            structured_replay(cert)
        with self.assertRaises(ValueError):
            compile_h(hamiltonian(2, F(4)), 2, 1)
        with self.assertRaises(ValueError):
            compile_h(hamiltonian(4, F(4)), 4, 5)


if __name__ == "__main__":
    unittest.main()
