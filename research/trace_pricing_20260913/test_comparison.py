"""Comparison must not interpolate unavailable resource/accuracy outcomes."""
from fractions import Fraction as F
import unittest

from research.trace_pricing_20260913.replay import best_within


class ResourceComparison(unittest.TestCase):
    def test_counts_and_elapsed_time_filter_completed_points(self):
        points = [{'case': 'seed', 'directions': 64, 'complete_at_seconds': 0., 'width_Ha': '7/1000'},
            {'case': 'small', 'directions': 80, 'complete_at_seconds': 30., 'width_Ha': '5/1000'},
            {'case': 'large', 'directions': 96, 'complete_at_seconds': 70., 'width_Ha': '2/1000'}]
        self.assertEqual(best_within(points, count=90)['case'], 'small')
        self.assertEqual(best_within(points, elapsed=69.999)['case'], 'small')
        self.assertEqual(best_within(points, elapsed=70.)['case'], 'large')
        self.assertEqual(best_within(points, count=80, elapsed=29.)['case'], 'seed')
        self.assertIsNone(best_within(points, count=63))

    def test_best_bound_not_merely_last_completed_trial(self):
        points = [{'case': 'better', 'directions': 80, 'complete_at_seconds': 30., 'width_Ha': '1/200'},
            {'case': 'later_weaker', 'directions': 96, 'complete_at_seconds': 60., 'width_Ha': '3/500'}]
        self.assertEqual(best_within(points, elapsed=100.)['case'], 'better')


if __name__ == '__main__':
    unittest.main()
