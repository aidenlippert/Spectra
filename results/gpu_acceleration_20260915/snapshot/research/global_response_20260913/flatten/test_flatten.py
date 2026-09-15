import unittest
from fractions import Fraction as F
from research.composable_response_20260913 import recursive
from .flatten import nested_action, flattened_action, schedule_counts, factored_response_constant


class FlattenTests(unittest.TestCase):
    def test_chebyshev_factored_response_identity(self):
        # F_k = 2p_k-d p_k^2 = c_k p_{2k}, checked at rational points.
        for k in (1,2,4):
            delta=F(2); M=F(10); c=(M+delta)/2; a=(M-delta)/2; z0=c/a
            ck=factored_response_constant(k,z0)
            def p(degree,d):
                from research.compact_response_20260913 import program
                z=(c-d)/a
                return (1-program.chebyshev(degree,z)/program.chebyshev(degree,z0))/d
            for d in (F(3),F(7),F(9)):
                lhs=2*p(k,d)-d*p(k,d)*p(k,d); rhs=ck*p(2*k,d)
                self.assertEqual(lhs,rhs)

    def setUp(self):
        self.H=[[F(3),F(1,5),F(1,7),F(1,11)],
                [F(1,5),F(4),F(1,3),F(1,13)],
                [F(1,7),F(1,3),F(5),F(1,17)],
                [F(1,11),F(1,13),F(1,17),F(6)]]
        self.first={'delta_Ha':'2','M_Ha':'6','order':4}
        self.second={'delta_Ha':'1','M_Ha':'8','order':3}
        self.P1=lambda v:[v[0],v[1],v[2],F(0)]
        self.Q1=lambda v:[F(0),F(0),F(0),v[3]]
        self.P2=lambda v:[v[0],v[1],F(0),F(0)]
        self.Q2=lambda v:[F(0),F(0),v[2],F(0)]
    def test_flatten_preserves_noncommuting_order(self):
        v=[F(1),F(2),F(3),F(0)]
        expected=nested_action(lambda x:[sum(a*b for a,b in zip(row,x)) for row in self.H],self.P1,self.Q1,self.P2,self.Q2,v,self.first,self.second)
        got,stats=flattened_action(lambda x:[sum(a*b for a,b in zip(row,x)) for row in self.H],self.P1,self.Q1,self.P2,self.Q2,v,self.first,self.second)
        self.assertEqual(got,expected)
        self.assertEqual(stats['unique_H_rhs'],stats['H_rhs'])
        self.assertEqual(schedule_counts(4,3)['full_K2'],63)
    def test_cse_shares_identical_rhs_only(self):
        # Zeroing the outer input makes several polynomial intermediates equal;
        # this exercises genuine sharing without assuming commutativity.
        zero=[F(0)]*4
        H=lambda x:[sum(a*b for a,b in zip(row,x)) for row in self.H]
        got,stats=flattened_action(H,self.P1,self.Q1,self.P2,self.Q2,zero,self.first,self.second)
        self.assertEqual(got,nested_action(H,self.P1,self.Q1,self.P2,self.Q2,zero,self.first,self.second))
        self.assertLess(stats['unique_H_rhs'],schedule_counts(4,3)['full_K2'])


if __name__=='__main__': unittest.main()
