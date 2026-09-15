"""Exact algebra, fixed-N, soundness, and resource-refusal regression tests."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product as choices
from pathlib import Path
from unittest.mock import patch
import json
import unittest

from experiments.marginal_symbolic import add, adj, mono, product, scale
from experiments.marginal_transfer_verify import apply_word
from research.collective_interference_20260913.collective import (
    check_lower, core_forms, digest, model, parameters, rayleigh, replay,
    spectator_energy, states)
from research.collective_interference_20260913.validation import full_polynomial
from research.collective_interference_20260913.scope import one_leg_rank, residual, check_residual_witness

ROOT=Path(__file__).resolve().parents[2]
CAMPAIGN=ROOT/'results/collective_interference_20260913/campaign'


def fixture(name='dispersed_L2'):
    path=CAMPAIGN/name
    return json.loads((path/'model.json').read_text()),json.loads((path/'certificate.json').read_text())


class CollectiveTests(unittest.TestCase):
    def test_collective_car_norm_is_multiplicity(self):
        bright=add(*(mono(((0,i),)) for i in range(7)))
        self.assertEqual(add(product(bright,adj(bright)),product(adj(bright),bright)),mono((),7))

    def test_spectator_minimum_preserves_total_number(self):
        p=parameters(model(2),2)
        for endpoint in ('lo','hi'):
            for n in range(9):
                options=[]
                for occupations in choices(*(range(g['count']) for g in p['groups'])):
                    if sum(occupations)+n!=p['particles']:continue
                    options.append(sum(g[endpoint]*(a if g['sign']>0 else g['count']-1-a)
                                       for g,a in zip(p['groups'],occupations)))
                self.assertEqual(spectator_energy(p,n,endpoint),min(options) if options else None)

    def test_exact_physical_embedding_with_bright_and_dark_modes(self):
        # Four bins, each with two physical modes. Retain all 70 N_core=4
        # states and occupy the two negative-bin dark modes (N_total=6).
        data=model(2);p=parameters(data,2)
        vectors=[{i:1} for i in range(4)];dark=[]
        for g in p['groups']:
            offset=4+(p['count'] if g['sign']>0 else 0)
            indices=[offset+j for j in range(g['index_first'],g['index_last']+1)]
            vectors.append({i:1 for i in indices})
            if g['sign']<0:dark.append({indices[0]:1,indices[1]:-1})
        columns=[]
        for s in states(8,4):
            current={0:1}
            operators=[v for i,v in enumerate(vectors) if s>>i&1]+dark
            for vector in reversed(operators):
                output={}
                for state,c in current.items():
                    for i,a in vector.items():
                        target=apply_word(((1,i),),state)
                        if target:
                            t,sign=target;output[t]=output.get(t,0)+c*a*sign
                current={s:c for s,c in output.items() if c}
            self.assertTrue(all(s.bit_count()==6 for s in current))
            columns.append(current)
        for endpoint in ('lo','hi'):
            h,_=full_polynomial(data,2,endpoint);h=add(h,mono((),-p['filled_energy']))
            forms,_=core_forms(p,endpoint);f=next(f for f in forms if f['particles']==4)
            applied=[]
            for column in columns:
                result={}
                for w,c in h.items():
                    for s,a in column.items():
                        target=apply_word(w,s)
                        if target:
                            t,sign=target;result[t]=result.get(t,F(0))+c*a*sign
                applied.append(result)
            for i,left in enumerate(columns):
                for j,right in enumerate(columns):
                    overlap=sum(a*right.get(s,0) for s,a in left.items())
                    self.assertEqual(overlap,4*f['metric'][i] if i==j else 0)
                    entry=sum(a*applied[j].get(s,F(0)) for s,a in left.items())
                    self.assertEqual(entry,4*f['matrix'][i][j])

    def test_false_lower_is_rejected(self):
        data,cert=fixture();cert['lower_shifted']='100'
        with self.assertRaisesRegex(ValueError,'PSD'):replay(data,cert)

    def test_mismatched_model_is_rejected(self):
        data,cert=fixture();data['gap']='5'
        with self.assertRaisesRegex(ValueError,'binding'):replay(data,cert)

    def test_zero_upper_is_rejected(self):
        data,cert=fixture();cert['upper']['amplitudes']=[0]*len(cert['upper']['amplitudes'])
        with self.assertRaisesRegex(ValueError,'Nonzero'):replay(data,cert)

    def test_fractional_coordinates_are_rejected(self):
        data,cert=fixture();cert['upper']['amplitudes'][0]=0.5
        with self.assertRaisesRegex(ValueError,'integer'):replay(data,cert)

    def test_infeasible_upper_sector_is_rejected(self):
        data,cert=fixture('exact_L1');cert['upper']['particles']=2
        with self.assertRaisesRegex(ValueError,'total N'):replay(data,cert)

    def test_false_lower_trial_sector_is_rejected(self):
        data,cert=fixture('exact_L1');cert['lower_trial']['particles']=2
        with self.assertRaisesRegex(ValueError,'total N'):replay(data,cert)

    def test_pattern_budget_is_enforced(self):
        data,cert=fixture();cert['bins']=3
        with self.assertRaisesRegex(ValueError,'budget'):replay(data,cert)

    def test_gap_and_multiplicity_domain_is_enforced(self):
        for change in ({'gap':'0'},{'spread':'-1'},{'bath_length':0},{'bath_length':True},
                       {'active_number_offset':5}):
            with self.subTest(change=change),self.assertRaises(ValueError):parameters({**model(),**change},1)

    def test_large_replay_uses_only_retained_states(self):
        data,cert=fixture('dispersed_L64')
        original_states=states
        def bounded(m,n):
            self.assertLessEqual(m,8)
            return original_states(m,n)
        with patch('research.collective_interference_20260913.collective.states',side_effect=bounded), \
             patch('research.collective_interference_20260913.validation.full_polynomial',side_effect=AssertionError('Oracle called')):
            receipt=replay(data,cert)
        self.assertEqual(receipt['physical_modes'],8196)
        self.assertEqual(receipt['dark_modes_eliminated'],8190)

    def test_full_oracle_refuses_large_input(self):
        with self.assertRaisesRegex(ValueError,'capped'):full_polynomial(model(64))

    def test_lower_check_rejects_zero_pivot_coupling(self):
        with self.assertRaisesRegex(ValueError,'zero PSD pivot'):
            check_lower([{'matrix':[[F(0),F(1)],[F(1),F(2)]],'metric':[1,1]}],F(0))

    def test_uniform_ceiling_encloses_accepted_interval(self):
        data,cert=fixture();receipt=replay(data,cert)
        self.assertGreaterEqual(F(receipt['uniform_width_ceiling']),F(receipt['width']))
        self.assertGreaterEqual(F(receipt['lower_envelope_trial_gap']),0)
        self.assertEqual(F(receipt['maximum_bin_width']),F(1,2))

    def test_collective_density_is_removed_exactly_by_scope_grammar(self):
        ns=[mono(((1,i),(0,i))) for i in range(8)]
        active=add(*ns[:4]);bath=add(*ns[4:])
        h=add(scale(product(bath,bath),F(2,3)),scale(product(active,bath),F(5,7)))
        delta,_,_=residual(h,8,(0,1,2,3))
        self.assertFalse(delta)

    def test_rank_witness_does_not_reject_collective_elimination(self):
        h,p=full_polynomial(model(2));rank=one_leg_rank(h,12)
        self.assertEqual(rank['rank'],12)
        data,cert=fixture();self.assertGreaterEqual(F(replay(data,cert)['width']),0)

    def test_residual_witness_requires_exact_fixed_number_matrix_element(self):
        h=add(mono(((1,0),(0,1)),F(3,7)),mono(((1,1),(0,0)),F(3,7)))
        witness={'bra':1,'ket':2,'matrix_element':'3/7'}
        self.assertEqual(check_residual_witness(h,4,1,witness),F(3,7))
        with self.assertRaisesRegex(ValueError,'failed'):
            check_residual_witness(h,4,1,{**witness,'matrix_element':'4/7'})
        with self.assertRaisesRegex(ValueError,'fixed-N'):
            check_residual_witness(h,4,2,witness)


if __name__=='__main__':unittest.main()
