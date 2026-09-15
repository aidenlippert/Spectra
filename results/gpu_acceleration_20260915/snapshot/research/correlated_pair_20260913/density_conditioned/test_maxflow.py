import unittest
raise unittest.SkipTest('Rejected preliminary PH/flow implementation; replaced by test_algebra. See REJECTED_PROTOTYPES.md.')
from fractions import Fraction as F
from maxflow import feasible
def test_resonant_four_mode():
    assert feasible([F(1),F(1),F(1),F(1)],[(0,1,F(1)),(2,3,F(1))])['feasible']
    assert not feasible([F(0),F(0),F(0),F(0)],[(0,1,F(1))])['feasible']
def test_exact_cut():
    r=feasible([F(1,2),F(1,2)],[(0,1,F(1))]);assert r['flow']==F(1) and r['feasible']
if __name__=='__main__':test_resonant_four_mode();test_exact_cut();print('2 passed')
