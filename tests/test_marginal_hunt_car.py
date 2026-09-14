"""Independent occupation-basis checks of CAR polynomial arithmetic."""
from fractions import Fraction as F
from itertools import product
import unittest
from experiments.marginal_hunt_car import add, adj, mono, mul, run
from experiments.marginal_hunt_witness import family_receipt, local_identity, verify_dqg_blocks


def act(poly, state):
    """Direct Fock action; does not use CAR normal ordering."""
    out = {}
    for word, coefficient in poly.items():
        bits, amplitude = state, coefficient
        for creation, mode in reversed(word):
            occupied = (bits >> mode) & 1
            if occupied == bool(creation):
                amplitude = 0
                break
            amplitude *= (-1) ** ((bits & ((1 << mode) - 1)).bit_count())
            bits ^= 1 << mode
        out[bits] = out.get(bits, 0) + amplitude
    return {state: value for state, value in out.items() if value}


class CARTests(unittest.TestCase):
    def test_four_letter_words_match_independent_fock_action(self):
        letters = [(c, m) for c in (0, 1) for m in range(3)]
        for word in product(letters, repeat=4):
            raw = mono(word)
            reduced = mul(raw, mono(()))
            for state in range(8):
                self.assertEqual(act(raw, state), act(reduced, state))

    def test_cubic_and_cubic_linear_antipairs(self):
        examples = [
            add(mono(((0, 0), (0, 1), (0, 2)), F(2, 3)), mono(((0, 0), (0, 2), (0, 3)), -1)),
            add(mono(((1, 0), (0, 1), (0, 2)), 2), mono(((1, 1), (0, 2), (0, 3)), -1)),
            add(mono(((1, 0), (0, 0), (0, 2)), F(1, 3)), mono(((0, 2),), 2)),
        ]
        for b in examples:
            ba = adj(b)
            raw = add(*(mono(w + v, c * d) for w, c in b.items() for v, d in ba.items()),
                      *(mono(w + v, c * d) for w, c in ba.items() for v, d in b.items()))
            reduced = add(mul(b, ba), mul(ba, b))
            self.assertTrue(all(len(w) <= 4 for w in reduced))
            for state in range(16):
                self.assertEqual(act(raw, state), act(reduced, state))

    def test_single_squares_retain_degree_six(self):
        for kind in ("T1_charge_-3", "T2_charge_-1"):
            self.assertTrue(run(kind)["single_square_has_degree6"])
            self.assertEqual(run(kind)["degree6_terms"], 0)

    def test_corrupted_identity_fails_independent_action(self):
        valid = local_identity()
        corrupt = add(valid, mono((), F(1, 1000)))
        self.assertNotEqual(act(valid, 0), act(corrupt, 0))

    def test_dqg_blocks(self):
        self.assertEqual(verify_dqg_blocks()["G_entries_checked"], 1296)

    def test_scalable_family_and_exact_six_mode_ground_state(self):
        for r in (1, 2, 5, 16, 128):
            self.assertEqual(F(family_receipt(r)["closed_gap"]), F(r, 4))
        energies = []
        for state in range(64):
            if state.bit_count() == 3:
                populations = [(state & 7).bit_count(), (state >> 3).bit_count()]
                energies.append(sum(k * (k - 1) // 2 for k in populations))
        self.assertEqual(min(energies), 1)
        for invalid in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                family_receipt(invalid)


if __name__ == "__main__":
    unittest.main()
