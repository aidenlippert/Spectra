import unittest

from experiments.marginal_quartic_scan import probe_block


class QuarticScanTests(unittest.TestCase):
    def test_density_dressed_has_quartic_words(self):
        block, transform = probe_block(6, 'density_dressed')
        self.assertGreaterEqual(len(block['words']), 4)
        self.assertEqual(transform.shape[0], len(block['words']))

    def test_quartet_has_creator_and_annihilator_sectors(self):
        block, transform = probe_block(6, 'quartet')
        self.assertEqual(transform.shape[1], 2 * 15)
        self.assertTrue(any(all(c == 0 for c, _ in w) for w in block['words']))
        self.assertTrue(any(all(c == 1 for c, _ in w) for w in block['words']))


if __name__ == '__main__':
    unittest.main()
