"""Guard inherited descriptors and exact acceptance/target decisions."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

from research.spin_enrichment_20260913.campaign import SOURCE, TARGET, choose_best, validate_source


class Selection(unittest.TestCase):
    def row(self, case, lower, directions=80, accepted=True, elapsed=10.):
        return {'case': case, 'lower_Ha': str(lower), 'directions': directions,
            'accepted': accepted, 'complete_at_seconds': elapsed}

    def test_rejected_late_and_weaker_results_cannot_replace_source(self):
        rows = [self.row('late', '-1/1000', elapsed=361.),
            self.row('rejected', '0', accepted=False), self.row('weaker', '-1/50')]
        self.assertIsNone(choose_best(rows, F(0), 360., F(-1, 100)))
        rows.append(self.row('valid', '-9/1000'))
        self.assertEqual(choose_best(rows, F(0), 360., F(-1, 100)), 'valid')

    def test_target_is_exact_and_prefers_fewer_sufficient_directions(self):
        rows = [self.row('just_wide', -TARGET-F(1, 10**15), 72),
            self.row('at_target', -TARGET, 80), self.row('stronger_larger', -TARGET/2, 96)]
        self.assertEqual(choose_best(rows, F(0), 360.), 'at_target')


class InheritedSubspace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from research.spin_completion_20260913.discovery import Model
        root = Path(__file__).resolve().parents[2]/'results/molecular_collective_20260913/campaign/h6'
        cls.model = Model(json.loads((root/'fixture.json').read_text()), json.loads((root/'rank_10/tail.json').read_text()))
        cls.cert = json.loads((SOURCE/'certificate.json').read_text()); cls.span = json.loads((SOURCE/'span.json').read_text())

    def test_exact_saved_directions_match_frame_and_certificate(self):
        validate_source(self.model, self.cert, self.span)
        self.assertEqual(len(self.span), 64)

    def test_changed_binding_frame_or_direction_is_refused(self):
        cert = deepcopy(self.cert); cert['tail_sha256'] = 'wrong'
        with self.assertRaisesRegex(ValueError, 'binding failed'):
            validate_source(self.model, cert, self.span)
        span = deepcopy(self.span); span[0]['vector'][0] += 1
        with self.assertRaisesRegex(ValueError, 'do not match'):
            validate_source(self.model, self.cert, span)
        cert = deepcopy(self.cert); cert['anti_blocks'][0]['generators'][0][-1] ^= 1
        with self.assertRaisesRegex(ValueError, 'do not match'):
            validate_source(self.model, cert, self.span)


if __name__ == '__main__':
    unittest.main()
