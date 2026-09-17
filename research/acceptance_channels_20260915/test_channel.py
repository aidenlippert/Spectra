"""Independent exact occupation-action oracle for the supplied CAR identity."""
from fractions import Fraction
import unittest
from research.acceptance_channels_20260915.channel import check, supplied_channel


def act(poly, state, labels):
    out = {}
    for word, coefficient in poly.items():
        target, sign = state, 1
        for creation, mode in reversed(word):
            i = labels.index(mode)
            occupied = (target >> i) & 1
            if occupied == creation:
                sign = 0
                break
            sign *= (-1) ** ((target & ((1 << i)-1)).bit_count())
            target ^= 1 << i
        if sign:
            out[target] = out.get(target, Fraction(0)) + sign * coefficient
    return {k: v for k, v in out.items() if v}


def compose(left, right, state, labels):
    out = {}
    for middle, value in act(right, state, labels).items():
        for target, coefficient in act(left, middle, labels).items():
            out[target] = out.get(target, Fraction(0)) + value * coefficient
    return {k: v for k, v in out.items() if v}


class ChannelTest(unittest.TestCase):
    def test_all_128_local_basis_columns_against_literal_operator_action(self):
        B, _, P = supplied_channel()
        dagger = {tuple((1-c, i) for c, i in reversed(w)): value for w, value in B.items()}
        labels = sorted({i for w in B for _, i in w})
        self.assertEqual(len(labels), 7)
        for state in range(1 << len(labels)):
            expected = compose(dagger, B, state, labels)
            for target, value in compose(B, dagger, state, labels).items():
                expected[target] = expected.get(target, Fraction(0)) + value
            expected = {k: v for k, v in expected.items() if v}
            self.assertEqual(act(P, state, labels), expected)
            norm = sum(v*v for v in act(B, state, labels).values())
            norm += sum(v*v for v in act(dagger, state, labels).values())
            self.assertEqual(expected.get(state, 0), norm)

    def test_no_missing_dual_or_energy_claim(self):
        receipt = check()
        self.assertEqual(receipt['spatial_support'], [1, 2, 3, 4])
        self.assertFalse(receipt['contained_in_any_reported_old_window'])
        self.assertFalse(receipt['reported_negative_dual_value_replayed'])
        self.assertFalse(receipt['new_energy_bound'])


if __name__ == '__main__':
    unittest.main()
