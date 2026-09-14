import unittest
from fractions import Fraction as F

from experiments.marginal_signed_atoms import residual, verify_block,propose,price_atoms,replay,construct
from tests.test_marginal_clique_gap import fixture
from experiments.marginal_symbolic import add,mono,adj,product,encode
import contextlib
import io
import json
from pathlib import Path
import tempfile
from unittest.mock import patch


class SignedAtomTests(unittest.TestCase):
    def test_complete_support_pricing_and_exact_export(self):
        import numpy as np
        from math import comb
        b=np.eye(6); b[2:,2:]=np.full((4,4),-.4); np.fill_diagonal(b,1.)
        diagnostic={}
        found,minimum,count=price_atoms(b,set(),vector_mode='rational',pair_pricing='all_supports',diagnostics=diagnostic)
        self.assertEqual(count,sum(comb(6,k) for k in (2,3,4)))
        self.assertAlmostEqual(minimum,-.2)
        self.assertTrue(any(key[0]==(2,3,4,5) for key in found))
        self.assertEqual(diagnostic['complete_support_scan']['checked_total'],count)
        a=[[F(1) for j in range(4)] for i in range(4)]
        with contextlib.redirect_stdout(io.StringIO()):
            item,bound,d=propose(a,[1]*4,seconds=5,round_cap=2,vector_mode='rational',pair_pricing='all_supports')
        verify_block(a,item,bound,True)
        self.assertEqual(d['history'][0]['pricing']['complete_support_scan']['checked_total'],11)
        with self.assertRaises(ValueError): price_atoms(b,set(),pair_pricing='all_supports')

    def test_native_basis_roundtrip_and_semantic_refusal(self):
        a=[[F(1)-(F(2) if i==j else 0) for j in range(4)] for i in range(4)]
        with contextlib.redirect_stdout(io.StringIO()):
            item, bound, d = propose(a,[1]*4,seconds=2,round_cap=1)
        self.assertIsNotNone(d['native_basis'])
        basis=d['native_basis']
        with contextlib.redirect_stdout(io.StringIO()):
            _, resumed, rd = propose(a,[1]*4,seconds=2,round_cap=1,
                                     initial_item=item,
                                     initial_dictionary=d['atom_dictionary'],
                                     initial_basis=basis)
        self.assertAlmostEqual(float(resumed),float(bound),places=8)
        self.assertTrue(rd['initial_basis_applied'])
        bad=dict(basis,row_keys=basis['row_keys'][:-1])
        with self.assertRaises(ValueError):
            propose(a,[1]*4,seconds=1,round_cap=1,initial_item=item,
                    initial_dictionary=d['atom_dictionary'],initial_basis=bad)

    def test_native_basis_remaps_interleaved_fill_and_pending_atom_rows(self):
        a=[[F(1) for _ in range(4)] for _ in range(4)]
        a[0][0]=a[1][1]=2;a[0][1]=a[1][0]=0
        quad=((0,1,2,3),(1,1,1,1));pending=((0,1,2),(1,-1,1))
        with patch('experiments.marginal_signed_atoms.price_atoms',side_effect=[([quad],-1.,1),([pending],-1.,1)]),contextlib.redirect_stdout(io.StringIO()):
            item,bound,d=propose(a,[1]*4,seconds=5,round_cap=2,pair_mode='active')
        self.assertAlmostEqual(float(bound),0,places=8)
        basis=d['native_basis'];self.assertIsNotNone(basis)
        # Pair01 was appended after old atom rows. Rebuilt pair columns and
        # reversed atom order require semantic, not positional, remapping.
        self.assertNotEqual(basis['column_keys'][-1],'pair:2,3')
        reordered=dict(basis)
        for keys,statuses in [('column_keys','column_status'),('row_keys','row_status')]:
            reordered[keys]=list(reversed(basis[keys]));reordered[statuses]=list(reversed(basis[statuses]))
        with contextlib.redirect_stdout(io.StringIO()):
            result,resumed,rd=propose(a,[1]*4,seconds=5,round_cap=1,initial_item=item,
                initial_dictionary=list(reversed(d['atom_dictionary'])),initial_basis=reordered,pair_mode='active')
        self.assertTrue(rd['initial_basis_applied']);self.assertAlmostEqual(float(resumed),0,places=8)
        self.assertEqual(len(rd['native_basis']['row_keys']),len(basis['row_keys'])+1)
        verify_block(a,result,resumed)

    def test_basis_allows_more_basic_row_slacks_than_columns_and_rejects_bad_statuses(self):
        from itertools import combinations,product
        a=[[F(int(i==j)) for j in range(4)] for i in range(4)]
        initial={'metric_weights':[1]*4,'atom_scale':1,'atoms':[]}
        dictionary=[{'indices':list(g),'signs':[1]+list(signs)} for k in (3,4)
                    for g in combinations(range(4),k) for signs in product((-1,1),repeat=k-1)]
        with contextlib.redirect_stdout(io.StringIO()):
            item,bound,d=propose(a,[1]*4,seconds=5,round_cap=1,initial_item=initial,initial_dictionary=dictionary,pair_mode='active')
            result,resumed,rd=propose(a,[1]*4,seconds=5,round_cap=1,initial_item=item,
                initial_dictionary=d['atom_dictionary'],initial_basis=d['native_basis'],pair_mode='active')
        self.assertEqual(resumed,1);self.assertTrue(rd['initial_basis_applied'])
        basis=d['native_basis']
        for values in ([True]*len(basis['column_status']),[99]*len(basis['column_status']),[0]*len(basis['column_status'])):
            bad=dict(basis,column_status=values)
            with self.assertRaises(ValueError):
                propose(a,[1]*4,seconds=5,round_cap=1,initial_item=item,
                        initial_dictionary=d['atom_dictionary'],initial_basis=bad,pair_mode='active')
        bad=dict(basis,metric_weights=[2]*4)
        with self.assertRaises(ValueError):
            propose(a,[1]*4,seconds=1,round_cap=1,initial_item=item,
                    initial_dictionary=d['atom_dictionary'],initial_basis=bad)

    def test_support_balancing_selects_pairs_despite_more_negative_quadruples(self):
        import numpy as np
        b=np.eye(8)-np.ones((8,8))*.55
        low={};balanced={}
        found,_,_=price_atoms(b,set(),batch=12,vector_mode='rational',pair_pricing='all_pairs',diagnostics=low)
        spread,_,_=price_atoms(b,set(),batch=12,vector_mode='rational',pair_pricing='all_pairs',diagnostics=balanced,batch_selection='by_support')
        self.assertEqual(low['selected_by_support'],{'2':0,'3':0,'4':12})
        self.assertEqual(balanced['selected_by_support'],{'2':4,'3':4,'4':4})
        self.assertEqual(len(spread),len(set(spread)))
        self.assertTrue(all(float(np.asarray(v,dtype=float)@b[np.ix_(g,g)]@np.asarray(v,dtype=float))<0 for g,v in spread))
        with self.assertRaises(ValueError): price_atoms(b,set(),batch_selection='unknown')

    def test_complete_pair_pricing_finds_pair_outside_sampled_supports(self):
        import numpy as np
        b=np.eye(200)*10
        for i,j in ((20,21),(40,41),(60,61)):
            b[i,i]=b[j,j]=1;b[i,j]=b[j,i]=2
        b[150,150]=b[151,151]=.02;b[150,151]=b[151,150]=.03
        for i,start in ((150,160),(151,170)):
            for j in range(start,start+6): b[i,j]=b[j,i]=.1
        sampled,_,_=price_atoms(b,set(),batch=100000,vector_mode='rational')
        diagnostic={}
        complete,_,_=price_atoms(b,set(),batch=100000,vector_mode='rational',pair_pricing='all_pairs',diagnostics=diagnostic)
        self.assertNotIn((150,151),[g for g,v in sampled])
        self.assertIn((150,151),[g for g,v in complete])
        self.assertEqual(len(complete),len(set(complete)))
        self.assertEqual(diagnostic['by_support']['2']['supports'],19900)
        self.assertEqual(diagnostic['by_support']['2']['negative_eigenspaces'],4)
        self.assertTrue(diagnostic['all_pair_supports_checked'])
        with self.assertRaises(ValueError): price_atoms(b,set(),pair_pricing='all_pairs')

    def test_dictionary_preserves_final_priced_cuts_and_inactive_fill_pairs(self):
        a=[[F(int(i==j)) for j in range(4)] for i in range(4)]
        initial={'metric_weights':[1]*4,'atom_scale':1,'atoms':[]}
        atom=((0,1,2,3),(1,-1,1,-1))
        with patch('experiments.marginal_signed_atoms.price_atoms',return_value=([atom],-1.,1)),contextlib.redirect_stdout(io.StringIO()):
            item,bound,d=propose(a,[1]*4,seconds=5,round_cap=1,initial_item=initial,pair_mode='active')
        self.assertEqual(bound,1);self.assertEqual(d['candidate_atoms'],0)
        self.assertEqual(len(d['atom_dictionary']),1)
        # A cut that had zero weight and was never loaded in the prior last
        # round must activate all six zero-H pairs before its row is rebuilt.
        with contextlib.redirect_stdout(io.StringIO()):
            resumed,resumed_bound,r=propose(a,[1]*4,seconds=5,round_cap=1,initial_item=item,
                initial_dictionary=d['atom_dictionary'],pair_mode='active')
        self.assertEqual(r['initial_pairs'],6);self.assertEqual(r['initial_atoms'],1)
        self.assertEqual(resumed_bound,1);verify_block(a,resumed,1)
        for bad in ([d['atom_dictionary'][0]]*2,[dict(d['atom_dictionary'][0],weight=1)]):
            with self.assertRaises(ValueError): propose(a,[1]*4,initial_item=item,initial_dictionary=bad)

    def test_physical_rational_certificate_has_distinct_version_gate(self):
        b=add(mono(((0,1),(0,0))),mono(((0,3),(0,2)),4))
        c={'kind':'spin_rational_atom_complement_v1','modes':6,'particles':2,
           'hamiltonian':encode(product(adj(b),b)),'retained_states':[48],'target_lower':'0',
           'blocks':[{'states':[3,12],'metric_weights':[1,1],'atom_scale':1,
                      'atoms':[{'indices':[0,1],'amplitudes':[1,4],'amplitude_scale':1,'weight':1}]}]}
        c['blocks'] += [{'states':[s],'metric_weights':[1],'atom_scale':1,'atoms':[]} for s in (6,9,18,24,33,36)]
        self.assertEqual(F(replay(c)['target_lower']),0)
        c['kind']='spin_signed_atom_complement_v1'
        with self.assertRaises(ValueError): replay(c)

    def test_rational_amplitudes_reconstruct_and_use_quadratic_node_coefficients(self):
        a=[[F(-1),F(4)],[F(4),F(14)]]
        atom={'indices':[0,1],'amplitudes':[1,4],'amplitude_scale':4,'weight':16}
        initial={'metric_weights':[1,1],'atom_scale':1,'atoms':[atom]}
        r,m=residual(a,[1,1],[atom],1,True)
        self.assertEqual(r,[[F(-2),F(0)],[F(0),F(-2)]])
        with self.assertRaises(ValueError): residual(a,[1,1],[atom],1)
        with contextlib.redirect_stdout(io.StringIO()):
            item,bound,d=propose(a,[1,1],seconds=5,round_cap=2,initial_item=initial,pair_mode='active',vector_mode='rational')
        self.assertAlmostEqual(float(bound),-2,places=8)
        verify_block(a,item,bound,True)
        for values in ([0,4],[-1,4],[True,4]):
            bad=dict(atom,amplitudes=values)
            with self.assertRaises(ValueError): residual(a,[1,1],[bad],1,True)

    def test_local_eigenvector_pricing_returns_nonflat_rational_directions(self):
        import numpy as np
        b=np.array([[1.,-.5],[-.5,0.]])
        found,minimum,_=price_atoms(b,set(),vector_mode='rational')
        self.assertLess(minimum,0);self.assertTrue(found)
        group,v=found[0];vf=np.array(v,dtype=float)
        self.assertEqual(group,(0,1));self.assertNotEqual(abs(v[0]),abs(v[1]))
        self.assertLess(float(vf@b@vf),0)

    def test_active_pairs_add_fill_and_preserve_atom_multiplier_rows(self):
        # A=11^T+(e0-e1)(e0-e1)^T. The opposite edge cancels H01.
        a=[[F(1) for j in range(4)] for i in range(4)]
        a[0][0]=a[1][1]=2;a[0][1]=a[1][0]=0
        atom=((0,1,2,3),(1,1,1,1))
        results=[]
        for mode in ('active','all'):
            with patch('experiments.marginal_signed_atoms.price_atoms',side_effect=[([atom],-1.,1),([],0.,1)]),contextlib.redirect_stdout(io.StringIO()):
                item,bound,d=propose(a,[1]*4,seconds=5,round_cap=2,pair_mode=mode)
            verify_block(a,item,bound);self.assertAlmostEqual(float(bound),0,places=8)
            self.assertEqual(d['active_pairs'],6)
            self.assertTrue(any(x['indices']==list(atom[0]) for x in item['atoms']))
            results.append(d)
        self.assertEqual(results[0]['initial_pairs'],5)
        self.assertEqual(results[1]['initial_pairs'],6)

    def test_saved_h6_all_pair_signed_certificate(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6/signed_atom_adaptive_gap'
        c=json.loads((root/'certificate.json').read_text());r=replay(c)
        self.assertEqual(r,json.loads((root/'receipt.json').read_text()))
        self.assertEqual(r['unique_action_states'],368)
        self.assertGreater(r['target_lower_float'],-6.59)
        self.assertLess(r['target_lower_float'],-6.58)
        self.assertGreater(r['blocks'][1]['packing_threshold_lower_float'],-6.502)
        self.assertTrue(all(b['atom_maximum_support']==4 for b in r['blocks']))

    def test_native_proposal_handles_negative_threshold_and_all_pairs(self):
        a=[[F(1)-(F(2) if i==j else 0) for j in range(4)] for i in range(4)]
        with contextlib.redirect_stdout(io.StringIO()): item,bound,diagnostic=propose(a,[1]*4,seconds=5,round_cap=2)
        self.assertAlmostEqual(float(bound),-2,places=8)
        self.assertEqual(diagnostic['all_pairs'],6)
        verify_block(a,item,bound)
        import numpy as np
        b=np.array([[1,-.5,-.5,0],[-.5,0,0,0],[-.5,0,0,0],[0,0,0,0]])
        found,minimum,_=price_atoms(b,set())
        self.assertLess(minimum,0)
        self.assertTrue(found)

    def test_complete_physical_constructor_and_missing_q_rejection(self):
        c=fixture();c['kind']='spin_clique_metric_probe_v1';c['target_lower']='0'
        for b in c['blocks']: b['metric_weights']=[1]*len(b['states'])
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';source.write_text(json.dumps(c))
            with contextlib.redirect_stdout(io.StringIO()): r=construct(source,root/'out',seconds=5,round_cap=2)
            self.assertEqual(F(r['target_lower']),0)
            exported=json.loads((root/'out/certificate.json').read_text())
            with contextlib.redirect_stdout(io.StringIO()):
                warm=construct(source,root/'warm',seconds=5,round_cap=2,atom_source=root/'out/certificate.json',dictionary_source=root/'out/search_state.json')
            self.assertEqual(F(warm['target_lower']),0)
            self.assertTrue(warm['requested_target_certified'])
            bad=dict(exported,retained_states=[]);(root/'bad.json').write_text(json.dumps(bad))
            with self.assertRaises(ValueError): construct(source,root/'rejected',atom_source=root/'bad.json')
            dictionary=json.loads((root/'out/search_state.json').read_text())
            for field,value in [('retained_states',[]),('vector_mode','rational'),('hamiltonian',[]),('blocks',None),('blocks',[{}])]:
                changed=dict(dictionary,**{field:value});(root/'bad_dictionary.json').write_text(json.dumps(changed))
                with self.assertRaises(ValueError):
                    construct(source,root/'rejected',atom_source=root/'out/certificate.json',dictionary_source=root/'bad_dictionary.json')
            exported['blocks'].pop()
            with self.assertRaises(ValueError): replay(exported)

    def test_initial_exact_atoms_survive_exhausted_discovery_budget(self):
        a=[[F(1)-(F(2) if i==j else 0) for j in range(4)] for i in range(4)]
        initial={'metric_weights':[1]*4,'atom_scale':3,
                 'atoms':[{'indices':[0,1,2,3],'signs':[1]*4,'weight':3}]}
        item,bound,d=propose(a,[1]*4,seconds=1e-9,initial_item=initial)
        self.assertEqual(bound,-2)
        self.assertEqual(d['initial_atoms'],1)
        verify_block(a,item,-2)
        extra={'indices':[0,1,2],'signs':[1,-1,1]}
        dictionary=[extra,{k:v for k,v in initial['atoms'][0].items() if k!='weight'}]
        item,bound,d=propose(a,[1]*4,seconds=1e-9,initial_item=initial,initial_dictionary=dictionary)
        self.assertEqual(bound,-2);verify_block(a,item,-2)
        self.assertEqual(len(item['atoms']),1)
        self.assertEqual(F(item['atoms'][0]['weight'],item['atom_scale']),1)
        self.assertEqual(d['initial_atoms'],2);self.assertEqual(d['initial_proof_atoms'],1)
        with self.assertRaises(ValueError): propose(a,[1]*4,initial_item=initial,initial_dictionary=[extra])

    def test_mixed_sign_exact_residual(self):
        a = [[3, 1, -1], [1, 3, 1], [-1, 1, 3]]
        item = {'metric_weights': [1, 1, 1], 'atom_scale': 10,
                'atoms': [{'indices': [0, 1, 2], 'signs': [1, -1, 1], 'weight': 1}]}
        r, _ = residual(a, item['metric_weights'], item['atoms'], 10)
        self.assertEqual(r[0][1], F(11, 10))
        self.assertEqual(verify_block(a, item, F(-1))['atom_count'], 1)

    def test_same_support_different_signs_reconstructs_exactly(self):
        # H is the sum of two opposite-sign atoms plus a diagonal PSD margin.
        s = (1, 1, -1); t = (1, -1, 1)
        w = [F(2, 3), F(1, 3), F(1)]
        a = [[F(s[i]*s[j] + t[i]*t[j]) / (w[i]*w[j]) for j in range(3)] for i in range(3)]
        for i in range(3): a[i][i] += F(2)
        item = {'metric_weights': [2, 1, 3], 'atom_scale': 1,
                'atoms': [{'indices': [0, 1, 2], 'signs': list(s), 'weight': 1},
                          {'indices': [0, 1, 2], 'signs': list(t), 'weight': 1}]}
        # metric congruence leaves 2 W^2, hence the DD threshold is exactly 2.
        self.assertEqual(verify_block(a, item, F(2))['packing_threshold_lower'], '2')

    def test_bad_atom_fields_rejected(self):
        a = [[2, 1, 1], [1, 2, 1], [1, 1, 2]]
        base = {'metric_weights': [1, 1, 1], 'atom_scale': 10}
        for atom in ({'indices': [0, 1, 2], 'signs': [-1, 1, 1], 'weight': 1},
                     {'indices': [0, 1, 1], 'signs': [1, 1, 1], 'weight': 1},
                     {'indices': [0, 1, 2], 'signs': [1, 1, 1], 'weight': 0}):
            with self.assertRaises(ValueError):
                residual(a, [1, 1, 1], [atom], 10)

    def test_overstated_target_rejected(self):
        a = [[2, 1, 1], [1, 2, 1], [1, 1, 2]]
        item = {'metric_weights': [1, 1, 1], 'atom_scale': 10, 'atoms': []}
        with self.assertRaises(ValueError): verify_block(a, item, F(1))

    def test_metric_is_squared_in_bound(self):
        a = [[1, 4], [4, 16]]
        item = {'metric_weights': [4, 1], 'atom_scale': 1, 'atoms': []}
        # Congruence uses W=[1,1/4], so the transformed matrix is all ones.
        self.assertEqual(verify_block(a, item, F(-1))['packing_threshold_lower'], '0')


if __name__ == '__main__':
    unittest.main()
