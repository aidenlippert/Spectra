import copy
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import unittest

from experiments.marginal_signed_clique import verify_signed,dual_matrix,repair,replay,separator_pairing
from experiments.marginal_clique_gap import metric_matrix
from experiments.marginal_spin_reduction import SpinZeroOracle
from tests.test_marginal_clique_gap import fixture


class SignedCliqueTests(unittest.TestCase):
    def test_nonclique_signed_square_escapes_balanced_family(self):
        a=[[F(0),F(1),F(1)],[F(1),F(0),F(0)],[F(1),F(0),F(0)]]
        item={'metric_weights':[1,1,1],'node_weights':[1,0,0],'edge_weights':[0,0]}
        verify_signed(a,item)
        r=separator_pairing(a,item,[0,1,2],[1,1,1])
        self.assertEqual(F(r['unnormalized_dual_pairing']),-1)
        self.assertEqual(F(r['normalized_dual_atom_pairing']),F(-1,3))
        with self.assertRaises(ValueError): separator_pairing(a,item,[0,1,2],[1,-1,-1])
        with self.assertRaises(ValueError): separator_pairing(a,item,[0,1,1],[1,1,1])

    def test_independent_dual_matrix_trace_and_both_edge_signs(self):
        v=[1,2,3,4];gamma=F(2,7);weights=[12,6,4,3]
        a=[[F(v[i]*v[j])+(gamma if i==j else 0) for j in range(4)] for i in range(4)]
        item={'metric_weights':weights,'node_weights':[1,0,0,0],
              'edge_weights':[1,1,0,0,0,0]}
        r=verify_signed(a,item);b,m=dual_matrix(a,item)
        transformed,_=metric_matrix(a,weights)
        self.assertEqual(F(r['packing_threshold_upper']),gamma)
        trace=sum(transformed[i][j]*b[j][i] for i in range(4) for j in range(4))
        normalization=sum(m[i]*b[i][i] for i in range(4))
        self.assertEqual(trace/normalization,gamma)
        for i,j in combinations(range(4),2):
            for sign in (-1,1): self.assertGreaterEqual(b[i][i]+b[j][j]+2*sign*b[i][j],0)
        for size in (3,4):
            for group in combinations(range(4),size):
                self.assertGreaterEqual(sum(b[i][j] for i in group for j in group),0)
        bad=copy.deepcopy(item);bad['edge_weights'][3]=1
        with self.assertRaises(ValueError): verify_signed(a,bad)

    def test_identity_repair_restores_cliques_without_breaking_pair_caps(self):
        a=[[F(1) for j in range(4)] for i in range(4)]
        item,r=repair(a,[1]*4,[1,0,0,0],[0]*6,scale=12)
        self.assertEqual(r['identity_repair_integer'],6)
        verify_signed(a,item)
        b,_=dual_matrix(a,item)
        self.assertGreaterEqual(sum(sum(row) for row in b),0)
        with self.assertRaises(ValueError): repair(a,[1]*4,[1],[0]*6)

    def test_subset_support_extends_by_zero_and_rejects_p_states(self):
        c=fixture();c.pop('blocks');c.pop('complement_lower')
        c.update(kind='spin_signed_clique_dual_v1',states=[3,6,9,12],
                 metric_weights=[1]*4,node_weights=[1,0,0,0],edge_weights=[1,1,0,0,0,0],target_lower='1/100')
        r=replay(c)
        self.assertEqual(F(r['packing_threshold_upper']),0)
        self.assertEqual(r['unique_action_states'],4)
        self.assertTrue(r['target_excluded'])
        for mutation in ('p','duplicate','order','metric'):
            bad=copy.deepcopy(c)
            if mutation=='p': bad['states'][-1]=48
            elif mutation=='duplicate': bad['states'][-1]=9
            elif mutation=='order': bad['states'].reverse()
            else: bad.pop('metric_weights')
            with self.assertRaises(ValueError): replay(bad)

    def test_h6_signed_residual_obstruction_and_old_cover_refusal(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6'
        c=json.loads((root/'clique_signed_residual_obstruction/certificate.json').read_text())
        r=replay(c)
        self.assertTrue(r['target_excluded'])
        self.assertLess(r['packing_threshold_upper_float'],-6.74)
        self.assertGreater(r['packing_threshold_upper_float'],-6.75)
        self.assertEqual(r['unique_action_states'],168)
        old=json.loads((root/'clique_metric_width4_dual/certificate.json').read_text())
        b=old['blocks'][1];o=SpinZeroOracle(old)
        a=[[o.action(s).get(t,F(0)) for t in b['states']] for s in b['states']]
        with self.assertRaises(ValueError): verify_signed(a,b)


if __name__=='__main__': unittest.main()
