"""Independent M4 controls for sparse-map PSD block discovery."""
from __future__ import annotations

import json
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_symbolic import decode, verify
from research.certificate_scaling.adaptive_block_discovery import run


FIXTURE = Path("results/marginal_reynolds/m4_d4_exact/certificate.json")


class AdaptiveBlockControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = json.loads(FIXTURE.read_text())
        cls.h = decode(cls.source["hamiltonian"], 4, 4)

    def test_full_export_replay_and_dual_identity(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            receipt = run(self.h, 4, 2, out, width=2, rounds=2,
                          seconds=30, full=True, symmetry=False)
            cert = json.loads((out / "certificate.json").read_text())
            history = json.loads((out / "history.json").read_text())
        self.assertEqual(verify(cert)["lower"], receipt["lower"])
        self.assertAlmostEqual(history[0]["dual_identity"], 1.0, places=6)
        self.assertEqual(len(history), 1)  # full mode intentionally one solve

    def test_adaptive_supports_grow_and_replay_best_checkpoint(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            receipt = run(self.h, 4, 2, out, width=2, rounds=2,
                          seconds=30, full=False, symmetry=False, batch=8)
            history = json.loads((out / "history.json").read_text())
            cert = json.loads((out / "certificate.json").read_text())
        self.assertGreaterEqual(len(history), 2)
        self.assertTrue(all(history[i + 1]["supports_solved"] >= history[i]["supports_solved"]
                            for i in range(len(history) - 1)))
        self.assertEqual(verify(cert)["lower"], receipt["lower"])

    def test_symmetry_quotient_has_same_initial_numeric_bound(self):
        values = []
        receipts = []
        for symmetry in (False, True):
            with tempfile.TemporaryDirectory() as td:
                out = Path(td)
                receipt = run(self.h, 4, 2, out, width=2, rounds=2,
                              seconds=30, full=True, symmetry=symmetry)
                history = json.loads((out / "history.json").read_text())
            values.append(history[0]["numeric_lower"])
            receipts.append(receipt)
        self.assertAlmostEqual(values[0], values[1], delta=1e-4)
        self.assertTrue(receipts[1]["alpha_charge_quotient"])
        self.assertLess(receipts[1]["coefficient_rows"], receipts[0]["coefficient_rows"])

    def test_raw_snapshot_reconstructs_reported_residual_and_exact_eta(self):
        """The saved solver arrays must be independently auditable, not just metadata."""
        import numpy as np
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            receipt = run(self.h, 4, 2, out, width=2, rounds=1, seconds=30,
                          full=True, solver="CLARABEL", solver_eps=1e-7,
                          solver_seconds=20, save_raw=True, condition="row-column")
            history = json.loads((out / "history.json").read_text())[0]
            meta = json.loads((out / "raw_round_0.json").read_text())
            raw = np.load(out / "raw_round_0.npz")
            self.assertEqual(meta["condition"], "row-column")
            self.assertEqual(raw["x"].ndim, 1)
            self.assertEqual(raw["rem"].shape, (len(meta["rows"]),))
            self.assertEqual(len(meta["supports"]), len([k for k in raw.files if k.startswith("gram_")]))
            weights = np.asarray(meta["weights"], dtype=float)
            rem_l1 = float(weights @ np.abs(raw["rem"]))
            self.assertAlmostEqual(rem_l1, history["solver_residual_l1"], places=7)
            from experiments.marginal_symbolic import multiplier_basis,number_shift,mono,product,word_product
            rows=[tuple(tuple(letter) for letter in w) for w in meta['rows']]
            lookup={w:i for i,w in enumerate(rows)}
            # Reconstruct the actual polynomial, independently of the recorded
            # diagnostic and the sparse Gram-map implementation.
            coeff=np.zeros(len(rows));free=[mono(())]+[product(number_shift(4,2),p) for p in multiplier_basis(4,meta['ideal_body'])]
            for x,p in zip(raw['x'],free):
                for w,c in p.items():
                    if w in lookup:coeff[lookup[w]]+=x*float(c)
            for bi,block in enumerate(meta['supports']):
                words=[tuple(tuple(letter) for letter in w) for w in block['words']];q=raw[f'gram_{bi}']
                for i,left in enumerate(words):
                    dagger=tuple((1-c,p) for c,p in reversed(left))
                    for j,right in enumerate(words):
                        for w,c in word_product(dagger,right):
                            if w in lookup:coeff[lookup[w]]+=q[i,j]*float(c)
            rhs=np.array([float(self.h.get(w,0)) for w in rows])
            self.assertAlmostEqual(float(weights@abs(rhs-coeff)),history['raw_coefficient_residual_l1'],places=10)
            self.assertAlmostEqual(float(weights@abs(coeff+raw['rem']-rhs)),history['equality_violation_l1'],places=10)
            self.assertAlmostEqual(history["equality_violation_l1"], 0.0, delta=1e-6)
            self.assertEqual(verify(json.loads((out / "certificate.json").read_text()))["lower"],
                             receipt["lower"])
            self.assertAlmostEqual(history["exact_residual_l1"],
                                   history["export_numeric_residual_l1"], places=7)

    def test_row_column_conditioning_preserves_initial_bound(self):
        values = []
        for condition in ("none", "row-column"):
            with tempfile.TemporaryDirectory() as td:
                out = Path(td)
                run(self.h, 4, 2, out, width=2, rounds=1, seconds=30,
                    full=True, solver="CLARABEL", solver_eps=1e-7,
                    solver_seconds=20, condition=condition)
                values.append(json.loads((out / "history.json").read_text())[0]["numeric_lower"])
        self.assertAlmostEqual(values[0], values[1], delta=2e-5)

    def test_number_frame_exact_sector_polynomial_identities(self):
        from experiments.marginal_symbolic import mono,product,add,scale,number_shift
        for n in (2,3):
            shift=number_shift(4,n)
            for i in range(4):
                annih=mono(((0,i),));create=mono(((1,i),))
                numbers=[mono(((1,j),(0,j))) for j in range(4)]
                minus=add(scale(annih,n-1),*(scale(product(p,annih),-1) for p in numbers),product(annih,shift))
                plus=add(scale(create,n),*(scale(product(create,p),-1) for j,p in enumerate(numbers) if j!=i),product(create,shift))
                self.assertEqual(minus,{})
                self.assertEqual(plus,{})

    def test_mixed_linear_quotient_primal_dual_and_exact_replay(self):
        values = []
        for formulation in ("primal", "dual"):
            with tempfile.TemporaryDirectory() as td:
                out = Path(td)
                receipt = run(self.h, 4, 2, out, width=2, rounds=1, seconds=30,
                              full=True, family="mixed", ideal_body=2,
                              quotient_linear=True, anchor_identity=True,
                              formulation=formulation, solver="CLARABEL",
                              solver_eps=1e-7, solver_seconds=20)
                history = json.loads((out / "history.json").read_text())[0]
                cert = json.loads((out / "certificate.json").read_text())
                values.append(history["numeric_lower"])
                self.assertAlmostEqual(history["dual_identity"], 1.0, places=6)
                self.assertEqual(verify(cert)["lower"], receipt["lower"])
                self.assertTrue(receipt["quotient_linear"])
        self.assertAlmostEqual(values[0], values[1], delta=2e-5)


if __name__ == "__main__":
    unittest.main()
