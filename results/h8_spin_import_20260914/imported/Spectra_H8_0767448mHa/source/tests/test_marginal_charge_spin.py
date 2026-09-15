from fractions import Fraction as F
from itertools import combinations
import copy
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_charge_spin import ChargeSpin, cardinality_lower, quadratic_slice_gate, build, replay
from experiments.marginal_charge_polynomial import add as padd, multiply as pmultiply, scale as pscale
from experiments.marginal_coherent_tree import CoherentCharge
from experiments.marginal_symbolic import mono, add, product, adj, scale, encode
from tests.test_marginal_coherent_tree import fixture, direct_rows


def fixed_sign_fixture():
    data = fixture()
    n = [mono(((1,i),(0,i))) for i in range(8)]
    hopping = add(*(mono(((1,s),(0,2+s)), F(-1,3)) for s in (0,1)))
    hopping = add(hopping,adj(hopping))
    pair = mono(((1,0),(1,1),(0,3),(0,2)),F(1,5))
    data['hamiltonian'] = encode(add(*(scale(product(n[i],n[i+1]),4) for i in range(0,8,2)),
                                    hopping, scale(product(hopping,add(n[4],n[5])),F(1,3)),pair,adj(pair)))
    return data


class ChargeSpinTests(unittest.TestCase):
    def test_projected_psd_matches_balanced_four_spin_extrema(self):
        poly={0:F(7)}
        z=[{0:F(-1),1<<i:F(2)} for i in range(4)]
        for i,j in combinations(range(4),2):
            poly=padd(poly,pscale(pmultiply(z[i],z[j]),F((-1)**i*(i+j+1),7)))
        values=[sum((v for m,v in poly.items() if m&mask==m),F(0))
                for pair in combinations(range(4),2) for mask in [sum(1<<i for i in pair)]]
        gate=quadratic_slice_gate(poly,4,min(values)-F(1,100))
        self.assertEqual(gate['matrix_dimension'],3)
        self.assertIsNone(quadratic_slice_gate(poly,4,min(values)+F(1,100)))
        self.assertIsNone(quadratic_slice_gate({1:F(1)},4,F(-10)))
        self.assertIsNone(quadratic_slice_gate({7:F(1)},4,F(-10)))
        self.assertEqual(quadratic_slice_gate({0:F(7)},6,F(6))['matrix_dimension'],5)
        with self.assertRaises(ValueError): quadratic_slice_gate({1<<4:F(1)},4,F(-10))

    def test_build_replay_and_reject_overclaims(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory); source=p/'source.json'; source.write_text(json.dumps(fixed_sign_fixture()))
            with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
                result=build(source,p/'proof','2')
                c=json.loads((p/'proof/certificate.json').read_text())
                self.assertEqual(replay(c),result)
            self.assertEqual(result['covered_Q_configurations'],30)
            self.assertEqual(result['charge_patterns'],18)
            for key,value in [('target_lower','100'),('target_lower',2),('max_charge_patterns',17)]:
                bad=copy.deepcopy(c);bad[key]=value
                with self.assertRaises(ValueError):replay(bad)

    def test_energy_and_transfer_bind_actual_hamiltonian_and_valence_space(self):
        from tests.test_marginal_valence_reference import fixture as small_fixture
        from experiments.marginal_valence_reference import build as build_energy
        from experiments.marginal_spin_temple import replay as replay_energy
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);source=p/'source.json';source.write_text(json.dumps(small_fixture()))
            with contextlib.redirect_stdout(io.StringIO()):build_energy(source,p/'proof',upper_steps=2)
            c=json.loads((p/'proof/certificate.json').read_text());expected=replay_energy(c)['width']
            red=c['spin_symmetric_certificate']
            red['complement_atoms']={'kind':'valence_charge_spin_v1','metric_rule':{'family':'pair_square','factors':['1']*4},'max_charge_patterns':10000}
            self.assertEqual(replay_energy(c)['width'],expected)
            transfer=copy.deepcopy(c);r=transfer['spin_symmetric_certificate']
            r['complement_reference']={'hamiltonian':r['hamiltonian'],'complement_lower':r['complement_lower'],'complement_atoms':r.pop('complement_atoms')}
            self.assertEqual(replay_energy(transfer)['width'],expected)
            for mutate in [lambda r:r['complement_atoms'].update(hamiltonian=r['hamiltonian']),
                           lambda r:r.update(complement_lower='100'),
                           lambda r:r.update(retained_states=[r['retained_states'][0]])]:
                bad=copy.deepcopy(c);mutate(bad['spin_symmetric_certificate'])
                with self.assertRaises(ValueError):replay_energy(bad)

    def test_exact_rows_and_metric_sharing_for_all_charge_patterns(self):
        data=fixed_sign_fixture(); oracle=ChargeSpin(data); rows=direct_rows(CoherentCharge(data))
        covered=set(); patterns=list(oracle.patterns()); self.assertEqual(len(patterns),18)
        for q in patterns:
            with patch('experiments.marginal_determinant_tree.DeterminantOracle.action',side_effect=AssertionError('No actions')):
                poly,cost=oracle.compile(q)
            singles=[i for i,x in enumerate(q) if x==0]
            values=[]
            for state,value in rows.items():
                if tuple(((state>>(2*i))&3).bit_count()-1 for i in range(4)) != q: continue
                mask=sum(1<<j for j,i in enumerate(singles) if state & (1<<(2*i)))
                actual=sum((v for m,v in poly.items() if m&mask==m),F(0))
                self.assertEqual(actual,value); values.append(value); covered.add(state)
            self.assertEqual(len(values),cost['spin_assignments'])
            self.assertEqual(cost['spin_endpoint_evaluations'],0)
            self.assertLessEqual(cardinality_lower(poly,len(singles),len(singles)//2),min(values))
        self.assertEqual(covered,set(rows))

    def test_refusals_and_cardinality_lower(self):
        oracle=ChargeSpin(fixed_sign_fixture())
        for q in [(0,0,0,0),(1,0,0,0),(True,-1,0,0),(1,-1), (2,-2,0,0)]:
            with self.assertRaises(ValueError): oracle.compile(q)
        with self.assertRaises(ValueError): list(oracle.patterns(17))
        changing=ChargeSpin(fixture()); rejected=0
        for q in changing.patterns():
            try: changing.compile(q)
            except ValueError as exc:
                self.assertIn('changes sign',str(exc)); rejected+=1
        self.assertGreater(rejected,0)
        poly={0:F(2),1:F(-3),2:F(4),3:F(-10),5:F(1)}
        for occupied in range(4):
            values=[sum((v for m,v in poly.items() if m&mask==m),F(0))
                    for choice in combinations(range(3),occupied)
                    for mask in [sum(1<<i for i in choice)]]
            self.assertLessEqual(cardinality_lower(poly,3,occupied),min(values))


if __name__=='__main__': unittest.main()
