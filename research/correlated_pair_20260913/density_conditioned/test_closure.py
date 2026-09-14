import unittest
raise unittest.SkipTest('Rejected preliminary CAR implementation; replaced by test_algebra. See REJECTED_PROTOTYPES.md.')
from fractions import Fraction as F
from closure import factor,closure,check_all_terms
def test_quartic_and_degree_six_retained():
    f=factor(((1,0),),(1,1),{((1,2),):F(1)},F(1),F(1),F(1),{((1,3),(1,4),(0,3)):F(1)})
    c=closure(f);r=check_all_terms(c);assert r['degree6_terms']>0
    try:check_all_terms(c,True)
    except ValueError:return
    assert False
if __name__=='__main__':test_quartic_and_degree_six_retained();print('1 passed')
