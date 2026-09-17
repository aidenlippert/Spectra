"""Replay and tamper checks of actual eight-site proposal data, no full state."""
import contextlib,copy,io,json,unittest
from pathlib import Path
from research.constructive_compression_20260916.filter_checker import check

DATA=Path(__file__).resolve().parents[2]/'results/constructive_compression_20260916/filter64_exact'


def run(proof):
    with contextlib.redirect_stdout(io.StringIO()):return check(proof)


class FilterIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.target=json.loads((DATA/'proof.json').read_text())
        cls.loose=json.loads((DATA/'validation_loose_proof.json').read_text())

    def test_same_model_loose_bound_accepts(self):
        r=run(self.loose)
        self.assertTrue(r['accepted']);self.assertTrue(r['all_declared_steps_replayed'])
        self.assertEqual(r['enumerated_configurations'],0)
        self.assertEqual(r['lower'],'-3961527731/400000000')

    def test_requested_precision_rejects_from_actual_residuals(self):
        r=run(self.target)
        self.assertFalse(r['accepted']);self.assertEqual(r['checked_steps'],3)
        self.assertGreater(r['history'][-1]['error_floor_float'],2)

    def test_rejects_changed_seed(self):
        p=copy.deepcopy(self.target);p['states'][0]['tensors'][0][0][3]+=1
        with self.assertRaisesRegex(ValueError,'trace seed'):run(p)

    def test_rejects_changed_hamiltonian_binding(self):
        p=copy.deepcopy(self.target);p['states'][1]['fixture_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'Fixture binding'):run(p)

    def test_rejects_invalid_charge_flow(self):
        p=copy.deepcopy(self.target);p['states'][1]['tensors'][0][0][1]^=1
        with self.assertRaises(ValueError):run(p)

    def test_rejects_wrong_sector(self):
        p=copy.deepcopy(self.target);p['states'][1]['spin_counts']=[3,5]
        with self.assertRaises(ValueError):run(p)
