import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest

from experiments.marginal_clique_gap import margins,verify_cliques,propose_block,replay,construct,replay_flat,row_ceiling,replay_metric,metric_matrix
from experiments.marginal_factor_width import replay_transfer
from experiments.marginal_symbolic import add,mono,adj,product,encode


def fixture():
    b=add(mono(((0,1),(0,0))),mono(((0,3),(0,2))),
          mono(((0,3),(0,0))),mono(((0,2),(0,1)),-1))
    c={'kind':'spin_clique_complement_v1','modes':6,'particles':2,
       'hamiltonian':encode(product(adj(b),b)),'retained_states':[48],
       'complement_lower':'0','blocks':[{'states':[3,6,9,12],'clique_scale':1,
       'clique_atoms':[{'indices':[0,1,2,3],'weight':1}]}]}
    c['blocks'] += [{'states':[s],'clique_scale':1,'clique_atoms':[]} for s in (18,24,33,36)]
    return c


class CliqueGapTests(unittest.TestCase):
    def test_weighted_squares_use_congruence_identity_metric(self):
        v=[1,2,3,4];weights=[12,6,4,3];gamma=F(2,7)
        a=[[F(v[i]*v[j])+(gamma if i==j else 0) for j in range(4)] for i in range(4)]
        item={'metric_weights':weights,'clique_scale':1,
              'clique_atoms':[{'indices':[0,1,2,3],'weight':1}]}
        self.assertEqual(verify_cliques(a,gamma,item)['margin'],0)
        with self.assertRaises(ValueError): verify_cliques(a,gamma+F(1,1000),item)
        recipe,bound,_=propose_block(a,4,metric_weights=weights)
        self.assertAlmostEqual(float(bound),float(gamma),places=10)
        verify_cliques(a,bound,recipe)
        # Pull transformed square back explicitly; the residual is gamma*I.
        transformed,metric=metric_matrix(a,weights)
        diag,edges=margins(transformed,item['clique_atoms'],1)
        self.assertTrue(all(x==0 for x in edges.values()))
        self.assertEqual([x/m for x,m in zip(diag,metric)],[gamma]*4)
        self.assertEqual([[F(1)*F(12,weights[i])*F(12,weights[j])
            +(diag[i]/metric[i] if i==j else 0) for j in range(4)] for i in range(4)],a)
        for bad in ([0,1,1,1],[1,True,1,1],[1,1,1]):
            with self.assertRaises(ValueError): metric_matrix(a,bad)

    def test_three_and_four_coordinate_gap_on_rank_one_control(self):
        a=[[F(1) for j in range(4)] for i in range(4)]
        tri,b3,_=propose_block(a,3)
        quad,b4,_=propose_block(a,4)
        self.assertAlmostEqual(float(b3),-.5,places=9)
        self.assertEqual(b4,0)
        verify_cliques(a,b3,tri)
        verify_cliques(a,b4,quad)
        self.assertEqual(replay(fixture())['complement_lower'],'0')

    def test_exact_square_reconstruction_with_mixed_signs_and_edge_remainder(self):
        signs=[1,-1,1,-1];weight=F(1,5);gamma=F(1,7)
        a=[[weight*signs[i]*signs[j] for j in range(4)] for i in range(4)]
        for i in range(4): a[i][i]+=gamma+F(i+1,19)
        for i in range(4):
            for j in range(i+1,4):
                extra=F(i+j+1,31)
                a[i][j]+=extra*signs[i]*signs[j];a[j][i]=a[i][j]
                a[i][i]+=extra;a[j][j]+=extra
        item={'clique_scale':5,'clique_atoms':[{'indices':[0,1,2,3],'weight':1}]}
        diag,edges=margins(a,item['clique_atoms'],5)
        rebuilt=[[weight*signs[i]*signs[j] for j in range(4)] for i in range(4)]
        for (i,j),value in edges.items():
            rebuilt[i][i]+=value;rebuilt[j][j]+=value
            rebuilt[i][j]+=value*signs[i]*signs[j];rebuilt[j][i]=rebuilt[i][j]
        for i in range(4): rebuilt[i][i]+=diag[i]
        self.assertEqual(rebuilt,a)
        self.assertEqual(verify_cliques(a,gamma,item)['margin'],F(1,19))

    def test_rejects_false_signs_edge_budgets_and_thresholds(self):
        a=[[F(1) for j in range(4)] for i in range(4)]
        item={'clique_scale':1,'clique_atoms':[{'indices':[0,1,2,3],'weight':1}]}
        with self.assertRaises(ValueError): verify_cliques(a,F(1,100),item)
        for atom in ({'indices':[0,1,2,3],'weight':2},
                     {'indices':[0,1,1],'weight':1},
                     {'indices':[0,1,2,4],'weight':1},
                     {'indices':[0,1,2],'weight':0}):
            with self.assertRaises(ValueError): margins(a,[atom],1)
        with self.assertRaises(ValueError): margins(a,item['clique_atoms']*2,1)
        bad=copy.deepcopy(a);bad[0][1]=bad[1][0]=-1
        with self.assertRaises(ValueError): verify_cliques(bad,0,item)
        bad=fixture();bad['blocks'].pop()
        with self.assertRaises(ValueError): replay(bad)
        bad=fixture();bad['blocks'][0]['states'][0]=48
        with self.assertRaises(ValueError): replay(bad)

    def test_construct_ignores_old_factors_and_replays_new_squares(self):
        c=fixture()
        for block in c['blocks']:
            block.pop('clique_atoms');block.pop('clique_scale')
            block['factor']='invalid unused source factor'
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json'
            source.write_text(json.dumps({'spin_symmetric_certificate':c}))
            with contextlib.redirect_stdout(io.StringIO()): construct(source,root/'out',4)
            exported=json.loads((root/'out/certificate.json').read_text())
            self.assertEqual(F(replay(exported)['complement_lower']),0)
            self.assertTrue(all('factor' not in b for b in exported['blocks']))
            with self.assertRaises(ValueError): construct(source,root/'out',4)

    def test_constructor_binds_metric_to_hamiltonian_and_state_order(self):
        c=fixture();metric=copy.deepcopy(c);metric['kind']='spin_clique_metric_probe_v1'
        for b in metric['blocks']: b['metric_weights']=[1]*len(b['states'])
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';path=root/'metric.json'
            source.write_text(json.dumps({'spin_symmetric_certificate':c}))
            path.write_text(json.dumps(metric))
            with contextlib.redirect_stdout(io.StringIO()): construct(source,root/'out',4,path,5)
            exported=json.loads((root/'out/certificate.json').read_text())
            self.assertEqual(F(replay(exported)['complement_lower']),0)
            for change in ('states','hamiltonian','missing'):
                bad=copy.deepcopy(metric)
                if change=='states': bad['blocks'][0]['states'].reverse()
                elif change=='hamiltonian': bad['hamiltonian']=[]
                else: bad['blocks'][0].pop('metric_weights')
                path.write_text(json.dumps(bad))
                with self.assertRaises(ValueError): construct(source,root/change,4,path,5)
                self.assertFalse((root/change).exists())

    def test_flat_obstruction_does_not_exclude_unequal_amplitude_two_support_square(self):
        b=add(mono(((0,1),(0,0))),mono(((0,3),(0,2)),4))
        c={'kind':'spin_flat_factor_obstruction_v1','modes':6,'particles':2,
           'hamiltonian':encode(product(adj(b),b)),'retained_states':[48],
           'maximum_support':4,'target_lower':'0','row_state':3}
        r=replay_flat(c)
        self.assertEqual(F(r['threshold_ceiling']),F(-1,3))
        self.assertEqual(r['unique_action_states'],1)
        # H is already a two-coordinate positive rank-one matrix [1,4][1,4]^T.
        # Equal-magnitude support-four obstruction says nothing against that.
        c['maximum_support']=5
        with self.assertRaises(ValueError): replay_flat(c)
        c['maximum_support']=4;c['row_state']=48
        with self.assertRaises(ValueError): replay_flat(c)

    def test_saved_h6_clique_certificate_and_flat_transfer_obstruction(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6'
        c=json.loads((root/'clique_width4_gap/certificate.json').read_text())
        r=replay(c)
        self.assertEqual(r['complement_lower'],json.loads((root/'clique_width4_gap/receipt.json').read_text())['complement_lower'])
        self.assertEqual(r['unique_action_states'],368)
        self.assertTrue(all(x['maximum_atom_support']<=4 for x in r['blocks']))
        self.assertLess(r['complement_lower_float'],-7.5)
        obstruction=json.loads((root/'flat_width4_obstruction/certificate.json').read_text())
        q=replay_flat(obstruction)
        self.assertLess(q['threshold_ceiling_float'],-6.64)
        self.assertEqual(q['unique_action_states'],1)
        c=json.loads((root/'flat_width4_transfer_obstruction/certificate.json').read_text())
        q=replay_transfer(c)
        self.assertEqual(q['equal_magnitude_support_excluded'],4)
        self.assertNotIn('factor_width_excluded',q)
        self.assertGreater(q['threshold_deficit_float'],.35)

    def test_diagonal_metric_changes_only_the_necessary_gate(self):
        a=[[F(1),F(4)],[F(4),F(16)]]
        self.assertEqual(row_ceiling(a,[1,1]),F(-1,3))
        self.assertEqual(row_ceiling(a,[4,1]),F(2,3))
        # Positive row ceiling is not a spectral lower bound: det(A)=0.
        self.assertEqual(a[0][0]*a[1][1]-a[0][1]**2,0)
        for weights in ([1,0],[1,-1],[1],[1,True]):
            with self.assertRaises(ValueError): row_ceiling(a,weights)
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6/clique_metric_probe'
        c=json.loads((root/'certificate.json').read_text())
        r=replay_metric(c)
        self.assertEqual(r['unique_action_states'],368)
        self.assertTrue(all(b['necessary_row_gate_passes'] for b in r['blocks']))
        self.assertTrue(all(F(b['rescaled_threshold_ceiling'])>F(r['target_lower']) for b in r['blocks']))


if __name__=='__main__': unittest.main()
