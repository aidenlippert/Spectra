import unittest
from fractions import Fraction as F

from experiments.marginal_connected_dimer_moments import (
    compile_moments, cumulants, from_cumulants,
)
from experiments.marginal_clifford_moments import CliffordMomentOracle
from experiments.marginal_hubbard_polynomial import build
from experiments.marginal_transfer_verify import apply_word


def physical_action(h, vector):
    out = {}
    for state, coefficient in vector.items():
        for word, amplitude in h.items():
            image = apply_word(word, state)
            if image:
                target, sign = image
                out[target] = out.get(target, F(0)) + coefficient * amplitude * sign
    return {state: value for state, value in out.items() if value}


def dimer_boundary(sites):
    vector = {0: 1}
    for dimer in range(sites // 2):
        block = {0x9 << (4 * dimer): 1, 0x6 << (4 * dimer): -1}
        vector = {state | local: coefficient * amplitude
                  for state, coefficient in vector.items()
                  for local, amplitude in block.items()}
    return vector


def direct_normalized_moments(sites, U, t, order):
    certificate = build(sites, U, t)
    # Use the shared symbolic Hamiltonian decoding and independently apply its
    # CAR words to the physical boundary.
    oracle = CliffordMomentOracle(certificate)
    vector = dimer_boundary(sites)
    norm = sum(value * value for value in vector.values())
    powers = [vector]
    for _ in range(order):
        powers.append(physical_action(oracle.base.h, powers[-1]))
    return [sum(a * powers[k].get(state, 0) for state, a in vector.items()) / F(norm)
            for k in range(order + 1)]


class ConnectedDimerMomentTests(unittest.TestCase):
    def test_recurrence_is_exact_inverse(self):
        moments = [F(1), F(2), F(7), F(19), F(83), F(311)]
        self.assertEqual(from_cumulants(cumulants(moments)), moments)

    def test_small_chains_match_independent_car_contraction(self):
        for sites, U, t in [(2, 4, 1), (4, F(3, 2), F(2, 3)), (6, 0, 1)]:
            compiled, receipt = compile_moments(sites, U, t, 8)
            self.assertEqual(compiled, direct_normalized_moments(sites, U, t, 8))
            self.assertEqual(receipt['maximum_cluster_sites'], min(sites, 10))

    def test_zero_hamiltonian_has_only_normalization_moment(self):
        for sites in (2, 6, 20, 64):
            for U in (0, 4):
                moments, receipt = compile_moments(sites, U, 0, 8)
                self.assertEqual(moments, [1] + [0] * 8)
                self.assertEqual(receipt['cluster_count'], min(sites // 2, 5))
        moments, receipt = compile_moments(1000000, 4, 1, 0)
        self.assertEqual(moments, [1])
        self.assertEqual(receipt['cluster_count'], 0)
        self.assertEqual(receipt['total_source_frame_states'], 0)

    def test_odd_dimer_count_keeps_physical_norm_external(self):
        oracle = CliffordMomentOracle(build(6, 4, 1))
        normalized, _ = oracle.normalized_dimer_moments(8)
        physical, receipt = oracle.dimer_moments(8)
        scale = 2 ** (6 // 2)
        self.assertEqual(receipt['boundary_norm'], scale)
        self.assertEqual(physical, [[[x * scale for x in row] for row in matrix]
                                    for matrix in normalized])
        with self.assertRaisesRegex(ValueError, 'even number of Hadamards'):
            oracle.to_frame(dimer_boundary(6))

    def test_twelve_site_direct_clifford_comparison(self):
        compiled, receipt = compile_moments(12, 4, 1, 8)
        oracle = CliffordMomentOracle(build(12, 4, 1))
        direct, _ = oracle.normalized_dimer_moments(8)
        self.assertEqual([matrix[0][0] for matrix in direct], compiled)
        self.assertEqual(receipt['maximum_cluster_sites'], 10)

    def test_large_chain_uses_fixed_cluster_work(self):
        small, small_receipt = compile_moments(64, 4, 1, 8)
        large, large_receipt = compile_moments(10**6, 4, 1, 8)
        self.assertEqual(small_receipt['total_source_frame_states'],
                         large_receipt['total_source_frame_states'])
        self.assertEqual(large_receipt['cluster_count'], 5)
        self.assertEqual(large_receipt['maximum_cluster_sites'], 10)
        self.assertNotEqual(small, large)

    def test_refusal_gates(self):
        for args in [(3, 4, 1), (True, 4, 1), (10**9 + 2, 4, 1), (2, -1, 1),
                     (2, 4, -1), (2, 4, 1, 9), (2, 4.0, 1),
                     (2, 4, 1, True), (2, 4, 1, 1.0)]:
            with self.assertRaises(ValueError):
                compile_moments(*args)


if __name__ == '__main__':
    unittest.main()
