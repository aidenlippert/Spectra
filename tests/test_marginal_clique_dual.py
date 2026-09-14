import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from experiments.marginal_clique_dual import verify_block,propose,replay,replay_star
from experiments.marginal_clique_gap import propose_block,candidates
from tests.test_marginal_clique_gap import fixture


class CliqueDualTests(unittest.TestCase):
    def test_local_cover_complete_neighborhood_and_edge_constraints(self):
        c=fixture();c.pop('blocks');c.pop('complement_lower')
        c.update(kind='spin_clique_star_obstruction_v1',target_lower='1/100',row_state=3,
                 metric_states=[3,6,9,12],metric_weights=[1]*4,covered_edges=[[3,6],[3,9]])
        r=replay_star(c)
        self.assertEqual(F(r['packing_threshold_upper']),0)
        self.assertEqual(r['unique_action_states'],4)
        self.assertEqual(r['balanced_triangles_checked'],3)
        self.assertEqual(r['balanced_quads_checked'],1)
        for mutation in ('cover','neighborhood','metric','threshold'):
            bad=copy.deepcopy(c)
            if mutation=='cover': bad['covered_edges'].pop()
            elif mutation=='neighborhood': bad['metric_states'].pop()
            elif mutation=='metric': bad.pop('metric_weights')
            else: bad['target_lower']='0'
            with self.assertRaises(ValueError): replay_star(bad)

    def test_saved_h6_star_equals_full_dual_with_fewer_actions(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6'
        c=json.loads((root/'clique_metric_star_obstruction/certificate.json').read_text())
        r=replay_star(c)
        full=json.loads((root/'clique_metric_width4_dual/certificate.json').read_text())
        q=replay(full)
        self.assertEqual(r['packing_threshold_upper'],q['packing_threshold_upper'])
        self.assertEqual(r['unique_action_states'],50)
        self.assertEqual(r['covered_edges'],41)
        self.assertLess(r['packing_threshold_upper_float'],-6.79)
        self.assertTrue(q['target_excluded'])

    def test_exact_dual_pairing_and_nonconstant_metric(self):
        v=[1,2,3,4];gamma=F(2,7);weights=[12,6,4,3]
        a=[[F(v[i]*v[j])+(gamma if i==j else 0) for j in range(4)] for i in range(4)]
        # y=(3,0,0,0), z=(2,2,2,0,0,0) gives normalization 3.
        item={'metric_weights':weights,'node_weights':[3,0,0,0],'edge_weights':[2,2,2,0,0,0]}
        r=verify_block(a,item)
        self.assertEqual(F(r['packing_threshold_upper']),gamma)
        self.assertEqual(r['cliques_checked'],5)
        bad=copy.deepcopy(item);bad['edge_weights'][0]=1
        with self.assertRaises(ValueError): verify_block(a,bad)
        bad=copy.deepcopy(item);bad['node_weights']=[0]*4
        with self.assertRaises(ValueError): verify_block(a,bad)
        recipe,diagnostic=propose(a,weights,seconds=5)
        self.assertGreaterEqual(F(verify_block(a,recipe)['packing_threshold_upper']),gamma)
        self.assertAlmostEqual(diagnostic['packing_threshold_upper_float'],float(gamma),places=6)

    def test_full_dictionary_and_complete_physical_coverage(self):
        a=[[F(1) for j in range(4)] for i in range(4)]
        a[0][1]=a[1][0]=F(1,10**20)
        self.assertGreater(len(candidates(a,4,min_edge=0)),len(candidates(a,4)))
        c=fixture();c['kind']='spin_clique_packing_dual_v1';c['target_lower']='1/100'
        c.pop('complement_lower')
        for b in c['blocks']:
            n=len(b['states']);b.pop('clique_atoms');b.pop('clique_scale')
            b.update(metric_weights=[1]*n,node_weights=[3]+[0]*(n-1),
                     edge_weights=[2,2,2,0,0,0] if n==4 else [])
        r=replay(c)
        self.assertEqual(F(r['packing_threshold_upper']),0)
        self.assertTrue(r['target_excluded'])
        c['blocks'].pop()
        with self.assertRaises(ValueError): replay(c)


if __name__=='__main__': unittest.main()
