import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest

from experiments.marginal_determinant_tree import DeterminantOracle, GapFailure
from experiments.marginal_h6_fixture import interleaved_state
from experiments.marginal_reference_obstruction import replay
from experiments.marginal_sparse_molecular import run_source
from experiments.marginal_symbolic import decode
from experiments.marginal_transfer_verify import apply_word

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/marginal_h6'


class H6TransferTests(unittest.TestCase):
    def test_interleaving_signs_and_independent_integral_checks(self):
        for alpha in range(16):
            for beta in range(16):
                word = tuple((1, 2*i) for i in range(4) if alpha & (1 << i))
                word += tuple((1, 2*i+1) for i in range(4) if beta & (1 << i))
                self.assertEqual(interleaved_state(alpha, beta, 4), apply_word(word, 0))
        fixture = json.loads((DATA / 'fixture.json').read_text())
        reference = json.loads((DATA / 'reference_upper.json').read_text())
        oracle = DeterminantOracle(fixture)
        self.assertEqual((oracle.modes, oracle.particles, fixture['sector_dimension']), (12, 6, 924))
        for word in oracle.h:
            for spin in (0, 1):
                self.assertEqual(sum(2*c-1 for c, i in word if i % 2 == spin), 0)
        self.assertEqual(oracle.action(63)[63], F(reference['rhf_rational_diagonal']))
        self.assertEqual(oracle.upper(reference['independent_upper']), F(reference['upper']))
        self.assertLess(abs(float(F(reference['upper'])) - reference['electronic_fci_numerical']), 1e-8)

    def test_exact_retained_obstruction_and_false_claims(self):
        c = json.loads((DATA / 'retained_obstruction/certificate.json').read_text())
        fixture = json.loads((DATA / 'fixture.json').read_text())
        self.assertEqual(decode(c['hamiltonian'], 12, 4), decode(fixture['hamiltonian'], 12, 4))
        r = replay(c)
        self.assertGreater(r['retained_error_lower_bound_float'], .029)
        self.assertEqual(r['retained_error_lower_bound'], json.loads((DATA / 'retained_obstruction/receipt.json').read_text())['retained_error_lower_bound'])
        for update in ({'retained_floor': '0'}, {'retained_floor': '-100'},
                       {'retained_states': [63, 63]},
                       {'independent_upper': {'states': [63], 'amplitudes': [0]}}):
            bad = copy.deepcopy(c)
            bad.update(update)
            with self.assertRaises(ValueError): replay(bad)

    def test_updated_upper_does_not_make_the_row_gate_pass(self):
        fixture = json.loads((DATA / 'fixture.json').read_text())
        p = json.loads((DATA / 'sparse_p32/failure.json').read_text())['retained_states']
        d = json.loads((DATA / 'complement_gate_after_upper.json').read_text())
        oracle = DeterminantOracle(fixture)
        with self.assertRaises(GapFailure) as caught:
            oracle.cover(p, F(d['gamma']))
        self.assertEqual(caught.exception.state, d['counterexample_state'])
        self.assertEqual(caught.exception.lower, F(d['row_lower']))
        self.assertNotIn(caught.exception.state, p)

    def test_bare_larger_hamiltonian_reproduces_the_failed_selection(self):
        fixture = json.loads((DATA / 'fixture.json').read_text())
        bare = {k: fixture[k] for k in ('hamiltonian', 'modes', 'particles')}
        expected = json.loads((DATA / 'sparse_p32/failure.json').read_text())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'input.json').write_text(json.dumps(bare))
            with contextlib.redirect_stdout(io.StringIO()):
                run_source(root / 'input.json', root / 'result', selector='residual')
            got = json.loads((root / 'result/failure.json').read_text())
            self.assertEqual(got['status'], 'not_accepted')
            self.assertEqual(got['retained_states'], expected['retained_states'])
            self.assertEqual(got['independent_upper'], expected['independent_upper'])
            self.assertFalse((root / 'result/certificate.json').exists())


if __name__ == '__main__':
    unittest.main()
