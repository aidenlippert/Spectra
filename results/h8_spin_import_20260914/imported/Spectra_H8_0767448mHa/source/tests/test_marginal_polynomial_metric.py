import copy
import contextlib
import io
import json
from pathlib import Path
import tempfile
from fractions import Fraction as F
from math import lcm
import unittest
from unittest.mock import patch
from experiments.marginal_polynomial_metric import JointPolynomial, replay, check_nonsingleton_separator, replay_separator, complete_number_ideals
from experiments.marginal_symbolic import mono,product,scale,add,encode
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import spin_states
from tests.test_marginal_charge_spin import fixed_sign_fixture
from tests.test_marginal_coherent_tree import fixture


def data():
    c=fixed_sign_fixture();terms=[{'powers':[0]*4,'coefficient':10}]
    for i in range(4):
        powers=[0]*4;powers[i]=2;terms.append({'powers':powers,'coefficient':i+1})
    c['polynomial_metric']={'denominator':10,'terms':terms};return c


class PolynomialMetricTests(unittest.TestCase):
    def test_zero_valence_metric_removes_projector_without_losing_rows(self):
        c=data();c['polynomial_metric']={'denominator':2,'terms':[{'powers':[2 if i==j else 0 for i in range(4)],'coefficient':1} for j in range(4)]}
        ring=JointPolynomial(c);physical=SpinZeroOracle(c)
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
            _,poly,cost=ring.compile(F(2))
        self.assertTrue(cost['projector_free_zero_metric'])
        charges=lambda s:tuple(((s>>(2*i))&3).bit_count()-1 for i in range(4))
        for state in spin_states(physical):
            row=physical.action(state);expected=(row.get(state,F(0))-2)*ring.weight(charges(state))
            expected-=sum((abs(v)*ring.weight(charges(t)) for t,v in row.items() if t!=state and 1 in charges(t)),F(0))
            actual=sum((F(v,cost['numerator_scale']) for m,v in poly.items() if m&state==m),F(0))
            self.assertEqual(actual,expected)

    def test_zero_metric_keeps_local_chain_degree_bounded_at_eight_sites(self):
        modes=16;sites=8;n=[mono(((1,i),(0,i))) for i in range(modes)]
        h=add(*(scale(product(n[2*i],n[2*i+1]),4) for i in range(sites)),
              *(mono(((1,2*i+s),(0,2*(i+1)+s)),F(-1,3)) for i in range(sites-1) for s in (0,1)),
              *(mono(((1,2*(i+1)+s),(0,2*i+s)),F(-1,3)) for i in range(sites-1) for s in (0,1)))
        c={'modes':modes,'particles':sites,'hamiltonian':encode(h),
           'polynomial_metric':{'denominator':2,'terms':[{'powers':[2 if i==j else 0 for i in range(sites)],'coefficient':1} for j in range(sites)]}}
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
            _,_,cost=JointPolynomial(c).compile(F(2))
            c['polynomial_metric']={'denominator':1,'terms':[{'powers':[0]*sites,'coefficient':1}]}
            _,_,explicit=JointPolynomial(c).compile(F(2))
        self.assertTrue(cost['projector_free_zero_metric'])
        self.assertLessEqual(cost['numerator_degree'],4)
        self.assertEqual(explicit['numerator_degree'],8)

    def test_number_ideal_completion_is_an_exact_polynomial_identity(self):
        ring=JointPolynomial(data());poly={0:F(2,7),1:F(-5,3),3:F(1,11),21:F(4,9)}
        # Mask21 exceeds the fixed alpha population and is identically zero.
        poly=ring.multiply(poly,{0:1});remainder,ideals=complete_number_ideals(ring,poly)
        reconstructed=dict(remainder)
        for index,ideal in enumerate(ideals):
            shift={0:-ring.target,**{1<<i:1 for i in range(index,ring.modes,2)}}
            reconstructed=ring.add(reconstructed,ring.multiply(shift,ideal))
        self.assertEqual(reconstructed,poly)
        self.assertTrue(all(all((mask&s).bit_count()==ring.target for s in ring.spin_masks) for mask in remainder))
        for state in spin_states(ring.o.oracle):
            before=sum(v for m,v in poly.items() if m&state==m)
            after=sum(v for m,v in remainder.items() if m&state==m)
            self.assertEqual(before,after)

    def test_exact_nonsingleton_separator_and_residual_gate(self):
        ring=JointPolynomial(data());functional={'states':[15],'weights':['1']}
        result=check_nonsingleton_separator(ring,{0:-1},1,functional)
        self.assertEqual(result['functional_value'],'-1')
        self.assertTrue(result['residual_l1_allowance_covered'])
        self.assertGreater(result['positive_indicators_checked'],0)
        for proof in [{'states':[15],'weights':['2']},
                      {'states':[15,15],'weights':['1/2','1/2']},
                      {'states':[1],'weights':['1']},
                      {'states':[165],'weights':['1']},
                      {'states':[15,51],'weights':['2','-1']}]:
            with self.assertRaises(ValueError):check_nonsingleton_separator(ring,{0:-1},1,proof)
        with self.assertRaises(ValueError):check_nonsingleton_separator(ring,{0:1},1,functional)

    def test_separator_reconstructs_the_actual_polynomial(self):
        c=data();c.update(kind='joint_polynomial_nonsingleton_separator_v1',hamiltonian=encode({}),target_lower='1',
            polynomial_metric={'denominator':1,'terms':[{'powers':[0]*4,'coefficient':1}]},
            functional={'states':[15],'weights':['1']})
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
            self.assertEqual(replay_separator(c)['functional_value'],'-1')
        c['target_lower']='-1'
        with self.assertRaises(ValueError):replay_separator(c)

    def test_energy_and_reference_transfer_bind_hamiltonian_threshold_and_full_valence(self):
        from tests.test_marginal_valence_reference import fixture as small_fixture
        from experiments.marginal_valence_reference import build
        from experiments.marginal_spin_temple import replay as replay_energy
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=root/'source.json';source.write_text(json.dumps(small_fixture()))
            with contextlib.redirect_stdout(io.StringIO()):build(source,root/'proof',upper_steps=2)
            c=json.loads((root/'proof/certificate.json').read_text());expected=replay_energy(c)['width']
            red=c['spin_symmetric_certificate'];gamma=F(red['complement_lower'])-F(1,100)
            hopping=abs(next(F(t['coefficient']) for t in red['hamiltonian'] if t['word']==[[1,0],[0,2]]))
            den=lcm(gamma.denominator,hopping.denominator)
            red['complement_lower']=str(gamma)
            red['complement_atoms']={'kind':'joint_polynomial_metric_gap_v1',
                'polynomial_metric':{'denominator':1,'terms':[{'powers':[0,0],'coefficient':1}]},
                'weight_proof':{'denominator':1,'bound':1,'positive_indicators':[],'charge_indicators':[],'number_multipliers':[[],[]]},
                'numerator_proof':{'denominator':den,'bound':int((4-gamma)*den),'positive_indicators':[],
                    'charge_indicators':[{'required':0,'occupied':0,'weight':int((4+2*hopping)*den)}],
                    'number_multipliers':[[{'mask':0,'coefficient':int(-hopping*den)}] for _ in range(2)]}}
            self.assertEqual(replay_energy(c)['width'],expected)
            transfer=copy.deepcopy(c);r=transfer['spin_symmetric_certificate']
            r['complement_reference']={'hamiltonian':r['hamiltonian'],'complement_lower':r['complement_lower'],'complement_atoms':r.pop('complement_atoms')}
            self.assertEqual(replay_energy(transfer)['width'],expected)
            for mutate in [lambda r:r['complement_atoms'].update(hamiltonian=r['hamiltonian']),
                           lambda r:r['complement_atoms'].update(target_lower=r['complement_lower']),
                           lambda r:r.update(complement_lower='100'),
                           lambda r:r.update(retained_states=[r['retained_states'][0]]),
                           lambda r:r.update(hamiltonian=encode({}))]:
                bad=copy.deepcopy(c);mutate(bad['spin_symmetric_certificate'])
                with self.assertRaises(ValueError):replay_energy(bad)

    def test_global_positivity_certificate_and_refusals(self):
        c=data();n=[mono(((1,i),(0,i))) for i in range(8)]
        c['hamiltonian']=encode(add(*(scale(product(n[i],n[i+1]),4) for i in range(0,8,2))))
        c['polynomial_metric']={'denominator':1,'terms':[{'powers':[0]*4,'coefficient':1}]}
        c.update(kind='joint_polynomial_metric_gap_v1',target_lower='2',
                 weight_proof={'denominator':1,'bound':1,'positive_indicators':[],'charge_indicators':[],'number_multipliers':[[],[]]},
                 numerator_proof={'denominator':1,'bound':2,'positive_indicators':[],
                                  'charge_indicators':[{'required':0,'occupied':0,'weight':4}],'number_multipliers':[[],[]]})
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
            result=replay(c)
        self.assertEqual(result['weight_positivity']['lower'],'1')
        self.assertEqual(result['numerator_positivity']['lower'],'2')
        for mutate in [lambda c:c.update(target_lower='5'),
                       lambda c:c['numerator_proof']['charge_indicators'][0].update(weight=-4),
                       lambda c:c['weight_proof'].update(bound=0),
                       lambda c:c['polynomial_metric']['terms'][0].update(coefficient=-1)]:
            bad=copy.deepcopy(c);mutate(bad)
            with self.assertRaises(ValueError):replay(bad)

    def test_joint_numerator_matches_every_physical_row_and_excludes_target_valence(self):
        c=data();o=JointPolynomial(c);physical=SpinZeroOracle(c);states=spin_states(physical)
        with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
            u,k,cost=o.compile(F(2))
        def evaluate(poly,s,den):return sum((F(v,den) for m,v in poly.items() if m&s==m),F(0))
        def charge(s):return tuple(((s>>(2*i))&3).bit_count()-1 for i in range(4))
        saw_valence_source=False;saw_excluded_target=False
        for s in states:
            weight=o.weight(charge(s));self.assertEqual(evaluate(u,s,cost['metric_scale']),weight)
            row=physical.action(s);expected=(row.get(s,F(0))-2)*weight
            for t,value in row.items():
                if t==s:continue
                if 1 not in charge(t):
                    if 1 in charge(s):saw_excluded_target=True
                    continue
                expected-=abs(value)*o.weight(charge(t))
            self.assertEqual(evaluate(k,s,cost['numerator_scale']),expected)
            if 1 not in charge(s):saw_valence_source|=expected<0
        self.assertTrue(saw_valence_source);self.assertTrue(saw_excluded_target)
        self.assertEqual(cost['determinant_actions'],0)
        self.assertEqual(cost['occupation_endpoint_evaluations'],0)

    def test_slice_reduction_and_strict_metric_or_sign_refusals(self):
        c=data();o=JointPolynomial(c)
        self.assertEqual(o.multiply({(1<<0)|(1<<2):1},{1<<4:1}),{})
        for field,value in [('denominator',0),('denominator',True),('terms',[{'powers':[3,0,0,0],'coefficient':1}])]:
            bad=copy.deepcopy(c);bad['polynomial_metric'][field]=value
            with self.assertRaises(ValueError):JointPolynomial(bad)
        bad=fixture();bad['polynomial_metric']=c['polynomial_metric']
        with self.assertRaises(ValueError):JointPolynomial(bad).compile(F(2))


if __name__=='__main__':unittest.main()
