import copy
import json
from pathlib import Path
import tempfile
import unittest
from fractions import Fraction as F

from experiments.marginal_scaled_residual import scaled_bound, verify_block, replay, construct, refine_scaling, MAX_SCALING
from experiments.marginal_signed_atoms import replay as replay_atoms
from experiments.marginal_symbolic import add, mono, adj, product, encode


def physical_fixture():
    b=add(mono(((0,1),(0,0))),mono(((0,3),(0,2)),4))
    c={'kind':'spin_scaled_residual_complement_v1','modes':6,'particles':2,
       'hamiltonian':encode(product(adj(b),b)),'retained_states':[48],'target_lower':'0',
       'blocks':[{'states':[3,12],'metric_weights':[1,1],'atom_scale':1,
                  'atoms':[],'residual_scaling':[4,1]}]}
    c['blocks'] += [{'states':[s],'metric_weights':[1],'atom_scale':1,'atoms':[],
                    'residual_scaling':[1]} for s in (6,9,18,24,33,36)]
    return c


class ScaledResidualTests(unittest.TestCase):
    def test_nonconstant_original_metric_and_mixed_sign_pairs_reconstruct_exactly(self):
        t=[2,3,5]; w=[F(2,3),F(1,3),F(1)]; metric=[x*x for x in w]; gamma=F(-2)
        r=[[F(0) for _ in t] for _ in t]
        for i,j,value in [(0,1,F(1,5)),(0,2,F(-2,7)),(1,2,F(3,11))]: r[i][j]=r[j][i]=value
        for i in range(3): r[i][i]=gamma*metric[i]+sum(abs(r[i][j])*F(t[j],t[i]) for j in range(3) if j!=i)
        # Independently expand the implicit pair factors and diagonal margin.
        rebuilt=[[gamma*metric[i] if i==j else F(0) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(i+1,3):
                vector={i:t[j],j:(1 if r[i][j]>0 else -1)*t[i]}; coefficient=abs(r[i][j])/F(t[i]*t[j])
                for p,vp in vector.items():
                    for q,vq in vector.items(): rebuilt[p][q]+=coefficient*vp*vq
        self.assertEqual(rebuilt,r)
        a=[[r[i][j]/(w[i]*w[j]) for j in range(3)] for i in range(3)]
        item={'metric_weights':[2,1,3],'atom_scale':1,'atoms':[],'residual_scaling':t}
        self.assertEqual(F(verify_block(a,item,gamma)['scaled_residual_lower']),gamma)
        self.assertEqual(scaled_bound(r,metric,[7*x for x in t]),gamma)
        self.assertEqual(verify_block(a,item,gamma)['implicit_pair_count'],3)
        with self.assertRaises(ValueError): verify_block(a,item,gamma+F(1,1000))

    def test_scaling_rejects_nonpositive_wrong_length_bool_and_unbounded_values(self):
        r=[[F(1),F(4)],[F(4),F(16)]]
        for bad in (None,[],[1],[0,1],[-4,1],[True,1],[F(4),1],[MAX_SCALING+1,1]):
            with self.assertRaises(ValueError): scaled_bound(r,[F(1),F(1)],bad)

    def test_exact_refinement_increases_full_bound_without_changing_original_metric(self):
        r=[[F(1),F(1)],[F(1),F(4)]];m=[F(1),F(1)]
        t,bound,history=refine_scaling(r,m,[1,1])
        self.assertGreater(bound,F(2,3));self.assertEqual(bound,scaled_bound(r,m,t))
        self.assertTrue(any(h['common_factor']>1 for h in history))
        self.assertTrue(all(0<x<=MAX_SCALING for x in t))
        self.assertEqual(refine_scaling(r,m,[1,1],0),([1,1],F(0),[]))
        for diagonal in (F(1,10),F(1,3)):
            refined=refine_scaling([[diagonal]],[F(1)],[1])
            self.assertEqual(refined,([1],diagonal,[]))
            self.assertIsInstance(refined[1],F)
        with self.assertRaises(ValueError): refine_scaling(r,m,[1,1],65)

    def test_complete_physical_replay_version_and_coverage_gates(self):
        c=physical_fixture();r=replay(c)
        self.assertEqual(F(r['target_lower']),0);self.assertEqual(r['blocks'][0]['implicit_pair_count'],1)
        self.assertTrue(all(b['explicit_atom_count']==0 for b in r['blocks']))
        with self.assertRaises(ValueError): replay_atoms(c)
        changed=copy.deepcopy(c);changed['blocks'][0]['residual_scaling']=[1,1]
        with self.assertRaises(ValueError): replay(changed)
        changed=copy.deepcopy(c);changed['blocks'].pop()
        with self.assertRaises(ValueError): replay(changed)
        changed=copy.deepcopy(c);changed['kind']='spin_rational_atom_complement_v1'
        with self.assertRaises(ValueError): replay(changed)
        with self.assertRaises(ValueError): replay_atoms(changed)

    def test_constructor_preserves_stronger_scaling_and_checks_proposals(self):
        c=physical_fixture();old=copy.deepcopy(c);old['kind']='spin_rational_atom_complement_v1';old['target_lower']='-3'
        for b in old['blocks']: del b['residual_scaling']
        proposal=[{'dimension':len(b['states']),'scaling_weights':b['residual_scaling']} for b in c['blocks']]
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';seed=root/'scaling.json'
            source.write_text(json.dumps(old));seed.write_text(json.dumps(proposal))
            result=construct(source,root/'out',seed,target='0')
            self.assertTrue(result['requested_target_certified']);self.assertEqual(F(result['target_lower']),0)
            proposal[0]['scaling_weights']=[1,1];seed.write_text(json.dumps(proposal))
            repeated=construct(root/'out/certificate.json',root/'again',seed,target='0',refine_steps=0)
            self.assertEqual(F(repeated['target_lower']),0)
            self.assertFalse(json.loads((root/'again/proposal.json').read_text())[0]['accepted_scaling'])
            proposal.pop();seed.write_text(json.dumps(proposal))
            with self.assertRaises(ValueError): construct(source,root/'bad',seed)

    def test_fresh_numerical_scaling_uses_same_exact_physical_gate(self):
        old=physical_fixture();old['kind']='spin_rational_atom_complement_v1';old['target_lower']='-3'
        for b in old['blocks']: del b['residual_scaling']
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';source.write_text(json.dumps(old))
            result=construct(source,root/'out',target='0')
            self.assertEqual(F(result['target_lower']),0)
            self.assertTrue(result['requested_target_certified'])


if __name__=='__main__': unittest.main()
