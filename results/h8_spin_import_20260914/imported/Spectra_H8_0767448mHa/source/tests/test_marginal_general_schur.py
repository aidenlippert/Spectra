import copy
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb
from pathlib import Path
import unittest
import numpy as np

from experiments.marginal_general_schur import (embedding,apply_columns,gram,prepare,lower,pivots,replay,fixture,norm_bound,annihilating_moments,resolvent_self_energy)
from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_schur_transfer import collective_matrix
from experiments.marginal_symbolic import add,adj,canonical,mono,scale,encode,decode,transform
from experiments.marginal_transfer_verify import apply_word

ROOT=Path(__file__).resolve().parents[1]


class GeneralSchurTests(unittest.TestCase):
    def test_signed_embedding_exactly_intertwines_fermionic_reference(self):
        z=embedding();h0=hopping_polynomial(10,F(1,5));hz=apply_columns(h0,z);j=collective_matrix()
        self.assertEqual(gram(z,z),[[F(comb(5,i)) if i==k else 0 for k in range(6)] for i in range(6)])
        self.assertTrue(any(a<0 for column in z for a in column.values()))
        for k,column in enumerate(hz):
            expected={s:a*j[i][k] for i in range(6) for s,a in z[i].items() if j[i][k]}
            self.assertEqual(column,expected)

    def test_connected_fixture_breaks_all_pair_charges_and_exchange(self):
        for interaction in (False,True):
            h=fixture(interaction=interaction)
            for i in range(5):
                self.assertTrue(any(sum((1 if c else -1) for c,m in w if m in (i,i+5)) for w in h))
            self.assertNotEqual(transform(h,[(i+5)%10 for i in range(10)]),h)
            neighbors={i:set() for i in range(10)}
            for w in h:
                if len(w)==2 and w[0][1]!=w[1][1]:
                    a,b=w[0][1],w[1][1];neighbors[a].add(b);neighbors[b].add(a)
            visited={0};todo=[0]
            while todo:
                for target in neighbors[todo.pop()]-visited:visited.add(target);todo.append(target)
            self.assertEqual(len(visited),10)

    def test_nonzero_internal_perturbation_is_retained(self):
        h0=hopping_polynomial(10,F(1,5));z=embedding();base=gram(z,apply_columns(h0,z))
        cycle=prepare(fixture());mixed=prepare(fixture(interaction=True))
        self.assertEqual(cycle['projected_h'],base)
        self.assertNotEqual(mixed['projected_h'],base)

    def test_paired_norm_covers_all_degree_four_transition_patterns(self):
        for count in (1,2):
            for creators in combinations(range(4),count):
                for annihilators in combinations(range(4),count):
                    if creators==annihilators:continue
                    p=mono(tuple((1,i) for i in creators)+tuple((0,i) for i in annihilators))
                    h=canonical(add(p,adj(p)));matrix=np.zeros((16,16),dtype=int)
                    for word,c in h.items():
                        for s in range(16):
                            target=apply_word(word,s)
                            if target:matrix[target[0],s]+=int(c)*target[1]
                    squared=matrix@matrix
                    self.assertTrue(np.array_equal(squared,np.diag(np.diag(squared))))
                    self.assertTrue(set(np.diag(squared)).issubset({0,1}))
                    self.assertEqual(norm_bound(h,'hermitian_pairs'),1)
        self.assertEqual(norm_bound(mono((),F(-3)),'hermitian_pairs'),3)

    def test_automatic_recurrence_and_resolvent_on_two_eigenvalues(self):
        h=mono(((1,0),(0,0)));columns=[{1:F(1),2:F(1)}]
        data=annihilating_moments(h,columns)
        self.assertEqual(data['annihilator'],[F(0),F(-1),F(1)])
        self.assertEqual(resolvent_self_energy(data['annihilator'],data['moments'],F(-1),size=1),[[F(3,2)]])
        with self.assertRaises(ValueError):resolvent_self_energy(data['annihilator'],data['moments'],F(0),size=1)
        with self.assertRaises(ValueError):annihilating_moments(h,columns,max_degree=1)
        self.assertEqual(annihilating_moments(h,[{}])['annihilator'],[F(1)])

    def test_resolvent_matches_full_reference_inverse(self):
        states=[s for s in range(1024) if s.bit_count()==5]
        def dense(columns):return np.array([[float(c.get(s,0)) for c in columns] for s in states])
        h0=hopping_polynomial(10,F(1,5));h=fixture(interaction=True);z=dense(embedding())
        matrix=dense(apply_columns(h0,[{s:F(1)} for s in states]))
        delta=dense(apply_columns(add(h,scale(h0,-1)),embedding()))
        w=delta-z@np.linalg.solve(z.T@z,z.T@delta)
        data=prepare(h,'hermitian_pairs',True);parameter=F(33,10)
        got=np.array(resolvent_self_energy(data['annihilator'],data['moments'],parameter),dtype=float)
        expected=w.T@np.linalg.solve(matrix-float(parameter)*np.eye(252),w)
        self.assertEqual(len(data['annihilator'])-1,6)
        self.assertTrue(np.allclose(got,expected,atol=1e-13,rtol=1e-10))

    def test_saved_exact_intervals_bind_actual_fixtures(self):
        for prefix,interaction in (('cycle',False),('mixed',True)):
            p=ROOT/'results/marginal_general_schur'/(prefix+'_1_1000_paired_resolvent')
            certificate=json.loads((p/'certificate.json').read_text())
            self.assertEqual(decode(certificate['hamiltonian'],10,4),fixture(F(1,1000),interaction))
            receipt=replay(certificate)
            self.assertGreater(receipt['width_float'],0)
            self.assertLess(receipt['width_float'],1e-6)
            self.assertEqual(receipt['resolvent_degree'],6 if interaction else 4)

    def test_rejects_false_lower_and_invalid_operator_classes(self):
        for h in (mono(((1,0),(0,1))),add(mono(((1,0),)),mono(((0,0),)))):
            with self.assertRaises(ValueError):prepare(h)
        with self.assertRaises(ValueError):prepare(fixture(),norm_method='unverified')
        with self.assertRaises(ValueError):prepare(fixture(),resolvent='false')
        data=prepare(fixture());self.assertIsNone(pivots(data,F(4)))

    def test_upper_witness_and_changed_hamiltonian_are_rechecked(self):
        p=ROOT/'results/marginal_general_schur/cycle_1_1000_paired_resolvent/certificate.json'
        certificate=json.loads(p.read_text());baseline=replay(certificate)
        changed=copy.deepcopy(certificate);changed['independent_upper']['numerical_energy']=-999
        self.assertEqual(replay(changed)['upper'],baseline['upper'])
        changed['independent_upper']['amplitudes']=[1]
        with self.assertRaises(ValueError):replay(changed)
        changed=copy.deepcopy(certificate)
        h=decode(changed['hamiltonian'],10,4);changed['hamiltonian']=encode(add(h,mono((),F(-1))))
        with self.assertRaises(ValueError):replay(changed)


if __name__=='__main__':unittest.main()
