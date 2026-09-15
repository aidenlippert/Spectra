import copy
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import unittest

from experiments.marginal_factor_width import replay, replay_transfer, replay_separator
from experiments.marginal_operator_norm import replay_q
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_temple import inertia
from experiments.marginal_symbolic import mono, add, product, adj, encode, scale


def fixture():
    b=add(mono(((0,1),(0,0))),mono(((0,3),(0,2))),
          mono(((0,3),(0,0))),mono(((0,2),(0,1)),-1))
    return {'kind':'spin_factor_width_obstruction_v1','modes':6,'particles':2,
        'hamiltonian':encode(product(adj(b),b)), 'retained_states':[48],
        'target_lower':'0','factor_width':3,
        'dual_weights':{'states':[3,6,9,12],'amplitudes':[1,1,1,1]}}


class FactorWidthTests(unittest.TestCase):
    def test_rank_one_positive_operator_requires_four_configuration_support(self):
        c=fixture()
        r=replay(c)
        self.assertEqual(F(r['normalized_dual_pairing']),F(-1,2))
        self.assertEqual(r['minimum_possible_maximum_support'],4)
        self.assertEqual(r['unique_action_states'],4)
        c['factor_width']=4
        with self.assertRaises(ValueError): replay(c)
        c['factor_width']=3
        c['target_lower']='-1/2'
        with self.assertRaises(ValueError): replay(c)

    def test_trace_matches_independent_dual_and_all_small_principal_minors(self):
        c=fixture()
        o=SpinZeroOracle(c)
        c['hamiltonian']=encode(add(scale(o.h,F(1,7)),
            mono(((1,0),(0,0)),F(1,100)),mono(((1,1),(0,1)),F(1,100))))
        c['dual_weights']['amplitudes']=[2,3,5,7]
        o=SpinZeroOracle(c)
        states=c['dual_weights']['states']; w=c['dual_weights']['amplitudes']
        a=[[o.action(t).get(s,F(0)) for t in states] for s in states]
        b=[[F(w[i]**2) if i==j else -F((a[i][j]>0)-(a[i][j]<0))*w[i]*w[j]/2
            for j in range(4)] for i in range(4)]
        self.assertTrue(any(x<0 for row in a for x in row))
        self.assertTrue(any(x>0 and x!=1 for row in a for x in row))
        for count in (1,2,3):
            for indices in combinations(range(4),count):
                self.assertEqual(inertia([[b[i][j] for j in indices] for i in indices])['negative'],0)
        trace=sum(a[i][j]*b[j][i] for i in range(4) for j in range(4))/sum(x*x for x in w)
        self.assertEqual(F(replay(c)['normalized_dual_pairing']),trace)
        self.assertLess(trace,0)

    def test_rejects_invalid_width_support_and_non_obstructing_weights(self):
        for update in ({'factor_width':True},{'factor_width':1},{'factor_width':65},
            {'dual_weights':{'states':[3,6,9,12],'amplitudes':[1,-1,1,1]}},
            {'dual_weights':{'states':[3,6,9,48],'amplitudes':[1,1,1,1]}},
            {'dual_weights':{'states':[3,3],'amplitudes':[1,1]}},
            {'dual_weights':{'states':[5],'amplitudes':[1]}},
            {'dual_weights':{'states':[3],'amplitudes':[1]}}):
            c=fixture(); c.update(update)
            with self.assertRaises(ValueError): replay(c)

    def test_saved_h6_obstruction_and_linear_threshold_ceiling(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6/factor_width3_obstruction'
        c=json.loads((root/'certificate.json').read_text())
        saved=json.loads((root/'receipt.json').read_text())
        r=replay(c)
        self.assertEqual(r['normalized_dual_pairing'],saved['normalized_dual_pairing'])
        self.assertEqual(r['witness_support'],96)
        self.assertEqual(r['unique_action_states'],96)
        self.assertLess(r['normalized_dual_pairing_float'],-.03)
        cutoff=F(r['target_lower'])+F(r['normalized_dual_pairing'])
        c['target_lower']=str(cutoff)
        with self.assertRaises(ValueError): replay(c)
        c['target_lower']=str(cutoff+F(1,10**15))
        self.assertEqual(F(replay(c)['normalized_dual_pairing']),-F(1,10**15))

    def test_saved_transfer_obstruction_binds_reference_norm_and_current_lower(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6/factor_width3_transfer_obstruction'
        c=json.loads((root/'certificate.json').read_text())
        r=replay_transfer(c)
        self.assertGreater(r['threshold_deficit_float'],.0019)
        self.assertEqual(F(r['exact_full_spin_sector_perturbation_norm']),F(1,25))
        self.assertEqual(r['current_interval_unique_source_states'],400)
        for change in ('sector','retained','norm'):
            bad=copy.deepcopy(c)
            if change=='sector': bad['current_interval']['spin_symmetric_certificate']['particles']=4
            if change=='retained': bad['current_interval']['spin_symmetric_certificate']['retained_states']=[]
            if change=='norm': bad['saturated_norm']['perturbation']=[]
            with self.assertRaises((ValueError,TypeError)): replay_transfer(bad)

    def test_compressed_norm_saturation_and_stronger_transfer_obstruction(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6'
        norm=json.loads((root/'hopping_q_norm_squares/certificate.json').read_text())
        r=replay_q(norm)
        self.assertEqual(F(r['operator_norm']),F(1,25))
        self.assertEqual(r['witness_support'],4)
        self.assertEqual(r['witness_action_states'],4)
        bad=copy.deepcopy(norm)
        bad['saturating_witness']['states']=[243,246,249,252]
        with self.assertRaises(ValueError): replay_q(bad)
        bad=copy.deepcopy(norm)
        bad['saturating_witness']['amplitudes']=[1,0,0,0]
        with self.assertRaises(ValueError): replay_q(bad)
        c=json.loads((root/'factor_width3_q_transfer_obstruction/certificate.json').read_text())
        r=replay_transfer(c)
        self.assertEqual(F(r['exact_q_compressed_perturbation_norm']),F(1,25))
        self.assertGreater(r['threshold_deficit_float'],.0019)
        c['saturated_norm']['retained_states']=[]
        with self.assertRaises(ValueError): replay_transfer(c)

    def test_four_coordinate_atom_separates_the_exact_dual(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_h6/factor_width4_separator'
        c=json.loads((root/'certificate.json').read_text())
        r=replay_separator(c)
        self.assertLess(F(r['normalized_dual_atom_pairing']),0)
        self.assertEqual(r['atom_support'],4)
        self.assertEqual(r['separator_check_source_states'],4)
        self.assertEqual(r['parent_obstruction_source_states'],96)
        for amplitudes in ([1,0,0,0],[1,1,1,1]):
            bad=copy.deepcopy(c)
            bad['separator']['amplitudes']=amplitudes
            with self.assertRaises(ValueError): replay_separator(bad)


if __name__=='__main__': unittest.main()
