import unittest
from experiments.v2_laws import CalibrationRecord, learn, predict, required_shots, verify_certificate

class LawTests(unittest.TestCase):
    def test_learn_and_predict(self):
        k = required_shots(3, 0.1, 0.01)
        recs = [CalibrationRecord((1,0,0), (1,)*k), CalibrationRecord((0,1,0), (0,)*k), CalibrationRecord((0,0,1), (0,)*k)]
        out = learn(recs, .1, .01)
        self.assertEqual(out.mask, (1,0,0)); self.assertTrue(out.certificate_valid)
        self.assertEqual(predict(out.mask, (1,1,1)), 1)
        self.assertTrue(verify_certificate(recs, out, .1, .01))

    def test_rank_deficient_abstains(self):
        recs = [CalibrationRecord((1,0), (1,)), CalibrationRecord((1,0), (1,))]
        out = learn(recs, 0.1, .1)
        self.assertIsNone(out.mask); self.assertFalse(out.certificate_valid); self.assertEqual(out.rank, 1)

    def test_validation(self):
        with self.assertRaises(ValueError): required_shots(2, .5, .1)
        with self.assertRaises(ValueError): learn([CalibrationRecord((1,), (0,0))], .1, .1)

    def test_overdetermined_and_inconsistent(self):
        recs = [CalibrationRecord((1,0), (1,)), CalibrationRecord((0,1), (0,)), CalibrationRecord((1,1), (1,))]
        self.assertIsNotNone(learn(recs, 0, .1).mask)
        bad = recs[:2] + [CalibrationRecord((1,1), (0,))]
        out = learn(bad, 0, .1); self.assertIsNone(out.mask); self.assertFalse(out.certificate_valid)

    def test_unseen_and_tamper(self):
        recs = [CalibrationRecord((1,0), (1,)), CalibrationRecord((0,1), (0,))]
        out = learn(recs, 0, .1)
        self.assertEqual(predict(out.mask, (1,1)), 1)
        tampered = recs[:1] + [CalibrationRecord((0,1), (1,))]
        self.assertFalse(verify_certificate(tampered, out, 0, .1))

if __name__ == '__main__': unittest.main()
