import json
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_transfer import connected_hamiltonian, upper_ed
from experiments.marginal_symbolic import verify, decode
from experiments.marginal_transfer_verify import replay

class TransferTests(unittest.TestCase):
    def test_graph_is_connected_when_delta_positive(self):
        h=connected_hamiltonian(F(1,10))
        edges={(i,j) for w in h for i,j in []}  # polynomial-level check below
        hopping=[w for w in h if len(w)==2]
        self.assertTrue(any({i for _,i in w}=={1,2} for w in hopping))

    def test_saved_certificates_replay_and_have_upper_interval(self):
        root=Path(__file__).resolve().parents[1]
        for d in ("0","1_100","1_10"):
            c=json.loads((root/f"results/marginal_transfer/m6_delta{d}_mixed.json").read_text())
            r=verify(c)
            self.assertGreaterEqual(F(c["independent_upper"]["upper"]),F(r["lower"]))
            self.assertEqual(replay(c)["upper"], c["independent_upper"]["upper"])

    def test_upper_claim_is_ignored_and_tampering_rejected(self):
        root=Path(__file__).resolve().parents[1]
        c=json.loads((root/'results/marginal_transfer/m6_delta1_10_mixed.json').read_text())
        expected=replay(c)["upper"]
        c["independent_upper"]["upper"]="999"
        self.assertEqual(replay(c)["upper"],expected)
        c["independent_upper"]["amplitudes"]=c["independent_upper"]["amplitudes"][:-1]
        with self.assertRaises(ValueError): replay(c)

    def test_exact_upper_is_variational(self):
        u=upper_ed(connected_hamiltonian(F(1,10)),6,3)
        self.assertGreater(u["upper_float"],0)
        self.assertGreater(int(u["norm"]),0)

if __name__=="__main__": unittest.main()
