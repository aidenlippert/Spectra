from fractions import Fraction as F
from research.global_response_20260913.interference.weighted_blocks import weighted_bound
def test_weighted_identity_matches_unweighted():
    n=[[F(1),F(2)],[F(3),F(4)]]; a=weighted_bound(n,(0,0,0,0))[0]
    assert a==max(map(sum,n))*max(map(sum,zip(*n)))
if __name__=='__main__': test_weighted_identity_matches_unweighted(); print('1 passed')
