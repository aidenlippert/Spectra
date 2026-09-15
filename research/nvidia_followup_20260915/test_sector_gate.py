"""The integer-spin argument must never silently extend to odd-electron cases."""
import unittest
from research.nvidia_followup_20260915.strict_replay import require_supported_sector


class SectorTests(unittest.TestCase):
    def test_even_campaign_cases(self):
        for m,n in ((8,4),(12,6),(14,10),(16,8),(20,10)):
            require_supported_sector({'modes':m,'particles':n})

    def test_odd_and_empty_spin_sectors_rejected(self):
        for m,n in ((6,3),(20,9),(7,4),(4,0),(4,4),(4,True),(4,2.0)):
            with self.assertRaises(ValueError):
                require_supported_sector({'modes':m,'particles':n})


if __name__=='__main__':
    unittest.main()
