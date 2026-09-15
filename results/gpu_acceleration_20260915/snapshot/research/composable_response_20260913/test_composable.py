from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
from unittest.mock import patch

from experiments.marginal_symbolic import decode
from research.response_consistency_20260913.separator import act
from research.compact_response_20260913 import program
from research.molecular_collective_20260913.core import matmul,transpose,eye
from research.composable_response_20260913 import coupling,spatial_gap,recursive,joint,finite_controls
from research.certificate_scaling.commutator_dual_witness import psd


class ComposableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data,cls.tail,_=program.load_case('h6')
        cls.response=json.loads((program.OUT/'h6_program.json').read_text())
        cls.cert=json.loads((recursive.OUT/'h6_recursive.json').read_text())
        cls.joint_cert=json.loads((recursive.OUT/'h6_joint_response.json').read_text())

    def test_joint_response_accepts_without_determinant_builder(self):
        with patch('research.compact_response_20260913.closure.blocks',side_effect=AssertionError('No enumeration allowed')):
            r=joint.check(self.data,self.tail,self.response,self.joint_cert)
        self.assertEqual(r['status'],'certified_response')
        self.assertGreater(F(r['second_sector_gap_Ha']),F(4,100))
        self.assertEqual(r['second_order'],119)
        self.assertLessEqual(F(r['second_residual_penalty_Ha']),F(1,10**6))
        self.assertEqual(r['many_body_states_enumerated'],0)
        self.assertFalse(r['terminal_retained_positivity_certified'])

    def test_joint_gap_refuses_invalid_exact_data(self):
        for key,value in (('occupation_multiplier','-1'),('chemical_potential',0.1),('lower_H_joint_Ha','100'),('last_spatial_orbitals',2.0)):
            c=deepcopy(self.joint_cert['joint_sector']);c[key]=value
            with self.assertRaises(ValueError):joint.check_gap(self.data,self.tail,c)
        c=deepcopy(self.joint_cert['joint_sector']);c['negative_part_majorant'][0][0]='-1'
        with self.assertRaises(ValueError):joint.check_gap(self.data,self.tail,c)

    def test_joint_partition_and_second_error_tampering_refused(self):
        for which in ('partition','penalty'):
            c=deepcopy(self.joint_cert)
            if which=='partition':c['joint_sector']['last_spatial_orbitals']=3
            else:c['second_response']['residual_penalty_Ha']='0'
            with self.assertRaises(ValueError):joint.check(self.data,self.tail,self.response,c)

    def test_nonpositive_gap_and_order_cap_are_not_success(self):
        for d in (F(0),F(-1,1000)):
            self.assertEqual(joint.outer_program({'delta':d,'M':F(10),'g':F(1)})['status'],'no_positive_second_gap')
        self.assertEqual(joint.outer_program({'delta':F(1,10000),'M':F(100),'g':F(100)})['status'],'order_cap_exhausted')

    def test_joint_lower_bound_survives_arbitrary_response(self):
        # H>=gamma I on the joint sector; a deliberately imperfect response
        # still gives K>=gamma I, carrying the induced correction intact.
        H=[[F(2),F(1,3),F(1,2)],[F(1,3),F(3),F(1,4)],[F(1,2),F(1,4),F(4)]]
        gamma=F(1);psd([[H[i][j]-(gamma if i==j else 0) for j in range(3)] for i in range(3)])
        B=[H[0][2],H[1][2]];D=H[2][2];X=[F(7,13),F(-2,5)]
        K=[[H[i][j]-B[i]*X[j]-X[i]*B[j]+D*X[i]*X[j] for j in range(2)] for i in range(2)]
        exact=[[H[i][j]-B[i]*B[j]/D for j in range(2)] for i in range(2)]
        psd([[K[i][j]-(gamma if i==j else 0) for j in range(2)] for i in range(2)])
        self.assertEqual([[K[i][j]-exact[i][j] for j in range(2)] for i in range(2)],
            [[D*(X[i]-B[i]/D)*(X[j]-B[j]/D) for j in range(2)] for i in range(2)])

    def test_matched_control_rejects_negative_witness_as_success(self):
        original_loads=json.loads
        def tampered(text,*a,**kw):
            c=original_loads(text,*a,**kw)
            if isinstance(c,dict) and c.get('kind')=='zero_response':
                c['blocks']+=c['failures'];c['failures']=[];c['all_blocks_pass']=True
            return c
        with patch('json.loads',side_effect=tampered):
            with self.assertRaisesRegex(ValueError,'classified correctly'):finite_controls.replay()

    def test_spatial_certificate_accepts_without_determinant_builder(self):
        with patch('research.compact_response_20260913.closure.blocks',side_effect=AssertionError('No enumeration allowed')):
            for cert in self.cert['sector_gaps']:
                r=spatial_gap.check(self.data,self.tail,cert)
                self.assertEqual(r['many_body_states_enumerated'],0)
                self.assertEqual(r['new_orbital_matrix_dimension'],5)
                self.assertEqual(r['largest_matrix_dimension'],21)

    def test_spatial_certificate_rejects_wrong_majorant_and_endpoint(self):
        for what in ('majorant','endpoint','binding'):
            cert=deepcopy(self.cert['sector_gaps'][0])
            if what=='majorant':cert['negative_part_majorant'][0][0]='-100'
            elif what=='endpoint':cert['lower_H_QQ_Ha']='100'
            else:cert['fixture_sha256']='wrong'
            with self.assertRaises(ValueError):spatial_gap.check(self.data,self.tail,cert)

    def test_local_blocks_match_independent_physical_CAR_actions(self):
        blocks,inputs,outputs=coupling.local_blocks(self.data);h=decode(self.data['hamiltonian'],12,4)
        for column,(local,rest) in enumerate(zip(inputs,(15,7,19))):
            source=rest+(local<<8); direct=act(h,{source:F(1)})
            expected={s:a for s,a in direct.items() if s>>(10)==3}
            actual={}
            for row,out in enumerate(outputs):
                for r,a in act(blocks[row][column],{rest:F(1)}).items():
                    actual[r+(out<<8)]=a
            self.assertEqual(actual,expected)

    def test_sqrt_is_an_outward_bound(self):
        for x in (F(0),F(1,7),F(2),F(9,16)):
            y=coupling.sqrt_up(x);self.assertGreaterEqual(y*y,x)
            self.assertLess((y-F(1,10**8))**2,x) if y else None

    def test_bounded_operator_expansion_refuses_expired_budget(self):
        with self.assertRaises(ValueError):coupling.bounded_product({():F(1)},{():F(1)},0,{'word_products':0})

    def test_recursive_counterexample_reproduces(self):
        r=recursive.check(self.data,self.tail,self.response,self.cert)
        self.assertFalse(r['second_positive_gap_certified'])
        self.assertLess(F(r['scalar_information_obstruction']['exact_retained_K_Ha']),0)
        self.assertEqual(r['many_body_states_enumerated'],0)

    def test_recursive_counterexample_and_sector_tampering_refused(self):
        for which in ('counter','sector'):
            c=deepcopy(self.cert)
            if which=='counter':c['scalar_information_obstruction']['retained_K_Ha']='0'
            else:c['sector_gaps'][1]['orbital']=5
            with self.assertRaises(ValueError):recursive.check(self.data,self.tail,self.response,c)

    def test_two_levels_equal_exact_congruences_and_charge_actions(self):
        H=[[F(3),F(1,5),F(1,7)],[F(1,5),F(4),F(1,3)],[F(1,7),F(1,3),F(5)]]
        response={'delta_Ha':'2','M_Ha':'6','order':4}
        calls=[0]
        def action(v):
            calls[0]+=1;return [sum(x*y for x,y in zip(row,v)) for row in H]
        P=lambda v:[v[0],v[1],F(0)];Q=lambda v:[F(0),F(0),v[2]]
        K=lambda v:recursive.retained_action(action,P,Q,v,response)
        p=recursive.scalar_polynomial(F(5),response)
        expected=[[H[i][j]-H[i][2]*H[2][j]*(2*p-5*p*p) for j in range(2)] for i in range(2)]
        first=[]
        for v in ([F(1),F(0),F(0)],[F(0),F(1),F(0)]):
            before=calls[0];first.append(K(v)[:2]);self.assertEqual(calls[0]-before,2*response['order']+1)
        self.assertEqual(transpose(first),expected)
        P2=lambda v:[v[0],F(0),F(0)];Q2=lambda v:[F(0),v[1],F(0)]
        second=recursive.retained_action(K,P2,Q2,[F(1),F(0),F(0)],response)
        p2=recursive.scalar_polynomial(expected[1][1],response)
        self.assertEqual(second[0],expected[0][0]-expected[0][1]**2*(2*p2-expected[1][1]*p2*p2))


if __name__=='__main__':unittest.main()
