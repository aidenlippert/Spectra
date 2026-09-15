from fractions import Fraction as F
from research.global_response_20260913.interference.weighted_blocks import weighted_bound

def test_exact_weighted_bound():
    n=[[F(1),F(2)],[F(3),F(4)]]
    v,r,c=weighted_bound(n,(-1,1,0,1))
    assert v==r*c and v>0

def test_refuse_nonpositive_or_wrong_dimension():
    try: weighted_bound([[F(1)]],(0,))
    except ValueError: pass
    else: assert False
    try: weighted_bound([[F(1)]],(-33,0))
    except ValueError: pass
    else: assert False

if __name__=='__main__': test_exact_weighted_bound(); test_refuse_nonpositive_or_wrong_dimension(); print('2 passed')
