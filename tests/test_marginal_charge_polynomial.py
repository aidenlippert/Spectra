import copy
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_charge_polynomial import envelope, replay, propose
from experiments.marginal_symbolic import add, mono, product, scale, encode


def diagonal_certificate(modes=4):
    h = add(*(scale(product(mono(((1, i), (0, i))),
                            mono(((1, i+1), (0, i+1)))), 4)
              for i in range(0, modes, 2)))
    return {'kind': 'valence_charge_polynomial_v1', 'modes': modes,
            'particles': modes//2, 'hamiltonian': encode(h), 'degree': 2,
            'b': '4', 'positive_indicators': [],
            'charge_indicators': [{'required': 0, 'occupied': 0, 'weight': '4'}],
            'number_multipliers': [[], []]}


class ChargePolynomialTests(unittest.TestCase):
    def test_large_sector_exact_charge_proof_never_calls_determinant_action(self):
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',
                   side_effect=AssertionError('Enumeration forbidden')):
            result = replay(diagonal_certificate(64))
        self.assertEqual(result['lower'], '4')
        self.assertEqual(result['residual_l1'], '0')
        self.assertEqual(result['determinant_actions'], 0)

    def test_termwise_envelope_is_operator_lower_even_with_interference(self):
        import numpy as np
        from experiments.marginal_transfer_verify import apply_word
        from experiments.marginal_symbolic import adj
        cert = diagonal_certificate()
        hop = mono(((1, 0), (0, 2)), F(-1, 3))
        density_hop = product(hop, mono(((1, 1), (0, 1))))
        cert['hamiltonian'] += encode(add(hop, adj(hop), scale(density_hop, -2),
                                          scale(adj(density_hop), -2)))
        oracle, f, _, _ = envelope(cert)
        matrix = np.zeros((16, 16))
        for s in range(16):
            for word, value in oracle.h.items():
                applied = apply_word(word, s)
                if applied:
                    t, sign = applied
                    matrix[t, s] += float(value*sign)
            matrix[s, s] -= float(sum(value for mask, value in f.items() if s & mask == mask))
        self.assertGreaterEqual(np.linalg.eigvalsh(matrix)[0], -1e-12)

    def test_lp_exports_exact_bound_and_preserves_previous_output(self):
        from experiments.marginal_symbolic import adj
        cert = diagonal_certificate()
        for spin in (0, 1):
            hop = mono(((1, spin), (0, 2+spin)), F(-1, 3))
            cert['hamiltonian'] += encode(add(hop, adj(hop)))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'source.json'
            source.write_text(json.dumps(cert))
            result = propose(source, root/'proof', degree=2)
            saved = json.loads((root/'proof/certificate.json').read_text())
            self.assertEqual(replay(saved), result)
            self.assertLessEqual(F(result['lower']), F(10, 3))
            self.assertGreater(F(result['lower']), F(10, 3)-F(1, 10**8))
            self.assertEqual(result['determinant_actions'], 0)
            with self.assertRaises(ValueError):
                propose(source, root/'proof', degree=2)

    def test_refuses_invalid_local_positivity_and_sectors(self):
        base = diagonal_certificate()
        mutations = [
            lambda c: c.update(particles=1),
            lambda c: c.update(degree=True),
            lambda c: c['charge_indicators'][0].update(weight='-1'),
            lambda c: c['charge_indicators'][0].update(required=1, occupied=2),
            lambda c: c['charge_indicators'][0].update(required=1, occupied=1),
            lambda c: c['number_multipliers'][0].append({'mask': 3, 'coefficient': '1'}),
            lambda c: c['hamiltonian'].extend(encode(mono(((1, 0), (0, 1))))),
        ]
        for mutation in mutations:
            c = copy.deepcopy(base); mutation(c)
            with self.assertRaises(ValueError): replay(c)

    def test_forged_bound_is_charged_to_exact_residual(self):
        cert = diagonal_certificate(); cert['b'] = '1000'
        self.assertEqual(replay(cert)['lower'], '4')


if __name__ == '__main__':
    unittest.main()
