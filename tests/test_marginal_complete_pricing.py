import unittest
import numpy as np

from experiments.marginal_complete_pricing import scan_supports


class CompletePricingTests(unittest.TestCase):
    def test_streaming_top_candidates_match_bruteforce_with_ties(self):
        from itertools import combinations
        for a in (np.full((8,8),-.4)+np.eye(8)*1.4,
                  np.diag(np.arange(8))-np.ones((8,8))):
            expected=[]
            for k in (2,3,4):
                for group in combinations(range(8),k):
                    value=float(np.linalg.eigvalsh(a[np.ix_(group,group)])[0])
                    if value < -1e-8: expected.append(value)
            result=scan_supports(a,batch_size=13,top_k=5)
            self.assertEqual(result['negative_support_count'],len(expected))
            np.testing.assert_allclose([x['eigenvalue'] for x in result['negative_supports']],sorted(expected)[:5],atol=1e-12)
            for k,stats in result['per_order_stats'].items():
                self.assertEqual(stats['checked'],result['possible_supports'][k])

    def test_counts_and_negative_direction(self):
        a = np.diag([-2.0, 1.0, 3.0, 4.0])
        result = scan_supports(a, batch_size=2, top_k=5)
        self.assertEqual(result["checked_supports"], {"2": 6, "3": 4, "4": 1})
        self.assertEqual(result["checked_total"], 11)
        self.assertAlmostEqual(result["minimum_eigenvalue"], -2.0)
        self.assertTrue(any(0 in x["support"] for x in result["negative_supports"]))

    def test_singular_leading_block_is_scanned(self):
        a = np.array([[0., 0., 0.], [0., 1., -2.], [0., -2., 1.]])
        result = scan_supports(a, max_support=3, tolerance=1e-9)
        self.assertLess(result["minimum_eigenvalue"], -0.9)
        self.assertGreaterEqual(result["near_boundary_count"], 1)

    def test_rejects_non_symmetric_and_nonfinite(self):
        with self.assertRaises(ValueError):
            scan_supports([[0., 1.], [0., 0.]])
        with self.assertRaises(ValueError):
            scan_supports(np.array([[np.nan, 0.], [0., 1.]]))

    def test_partial_final_batch_is_evaluated(self):
        a = np.zeros((6, 6))
        a[4, 4] = a[5, 5] = -1.0
        result = scan_supports(a, max_support=2, batch_size=7, top_k=4)
        self.assertEqual(result["checked_total"], 15)
        self.assertTrue(any(x["support"] == [4, 5] for x in result["negative_supports"]))


if __name__ == "__main__":
    unittest.main()
