from fractions import Fraction as F
from research.global_response_20260913.interference.frame_bound import frame_lower, certify

def test_scalar_and_radius():
    r=frame_lower([{'word':[],'coefficient':'2'}, {'word':[[1,0],[0,1]],'coefficient':'-3/2'}])
    assert r['lower']==F(1,2)

def test_refuses_target_beyond_bound():
    try: certify([{'word':[],'coefficient':'0'},{'word':[[1,0],[0,0]],'coefficient':'-1'}], '0')
    except ValueError: return
    assert False

if __name__=='__main__':
    test_scalar_and_radius(); test_refuses_target_beyond_bound(); print('2 passed')
