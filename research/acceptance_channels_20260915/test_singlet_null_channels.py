import unittest
from research.acceptance_channels_20260915.singlet_null_channels import channels
from research.acceptance_channels_20260915.test_channel import act


class SingletNullChannelsTest(unittest.TestCase):
    def test_exact_channels_on_complete_two_site_two_electron_singlet_basis(self):
        singlets=({3:1},{12:1},{9:1,6:-1})
        for creation in (False,True):
            for operator in channels(2,creation):
                for state in singlets:
                    result={}
                    for occupation,coefficient in state.items():
                        for target,value in act(operator,occupation,list(range(4))).items():
                            result[target]=result.get(target,0)+coefficient*value
                    self.assertFalse({k:v for k,v in result.items() if v})
        # The identity is sector-specific, not an assertion that the operator is zero.
        self.assertTrue(act(channels(2,False)[0],10,list(range(4))))


if __name__=='__main__':unittest.main()
