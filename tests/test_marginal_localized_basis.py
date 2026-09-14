import copy
import json
from pathlib import Path
import tempfile
import unittest
from fractions import Fraction as F

from experiments.marginal_localized_basis import rational_orthogonal, matrix_product, replay_rotation, sector_diagnostic, construct
from experiments.marginal_orbital_rotation import rotate
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_spin_constructor import select_reference
from experiments.marginal_symbolic import add, adj, mono, product, scale, encode


def fixture():
    number=[mono(((1,i),(0,i))) for i in range(4)]
    h=add(number[0],number[1],scale(number[2],-2),scale(number[3],-2))
    h=add(h,scale(product(number[0],number[1]),3),scale(product(number[2],number[3]),5))
    for spin in range(2):
        hop=mono(((1,spin),(0,2+spin)),F(1,7));h=add(h,hop,adj(hop))
    return {'modes':4,'particles':2,'hamiltonian':encode(h)},h


class LocalizedBasisTests(unittest.TestCase):
    def test_rational_givens_handles_rotations_reflections_and_pi_chart(self):
        cases=[[[F(3,5),F(-4,5)],[F(4,5),F(3,5)]],[[F(0),F(1)],[F(1),F(0)]],
               [[F(-1),F(0)],[F(0),F(-1)]]]
        for target in cases:
            u,d=rational_orthogonal(target)
            self.assertEqual(u,target)
            self.assertEqual(matrix_product(u,[list(r) for r in zip(*u)]),[[1,0],[0,1]])

    def test_random_dense_target_is_approximated_with_exact_orthogonality(self):
        import numpy as np
        q,_=np.linalg.qr(np.random.default_rng(73).normal(size=(5,5)))
        u,d=rational_orthogonal(q)
        self.assertLess(d['maximum_target_entry_error'],1e-4)
        self.assertEqual(matrix_product(u,[list(r) for r in zip(*u)]),[[int(i==j) for j in range(5)] for i in range(5)])
        for bad in ([[1,1],[0,1]],[[float('nan')]],[[1,0,0]]):
            with self.assertRaises(ValueError): rational_orthogonal(bad)
        with self.assertRaises(ValueError): rational_orthogonal([[1]],denominator=True)

    def test_exact_forward_inverse_spin_and_tampering_gates(self):
        data,h=fixture();u=[[F(3,5),F(4,5)],[F(-4,5),F(3,5)]]
        c=dict(data,kind='rational_spatial_rotation_v1',rotation=[[str(x) for x in r] for r in u],
               rotated_hamiltonian=encode(rotate(h,u)))
        r=replay_rotation(c);self.assertTrue(r['inverse_exact']);self.assertTrue(r['spin_symmetry_valid'])
        bad=copy.deepcopy(c);bad['rotation'][0][0]='1'
        with self.assertRaises(ValueError): replay_rotation(bad)
        bad=copy.deepcopy(c);bad['rotated_hamiltonian']=data['hamiltonian']
        with self.assertRaises(ValueError): replay_rotation(bad)

    def test_reference_seed_is_physical_and_selected_independently(self):
        data,h=fixture();oracle=SpinZeroOracle(data)
        default,_,_=select_reference(oracle,2)
        chosen,_,_=select_reference(oracle,2,seed=12)
        self.assertEqual(default[0],3);self.assertEqual(chosen[0],12)
        for bad in (True,1,5):
            with self.assertRaises(ValueError): select_reference(oracle,2,seed=bad)
        a,b=sector_diagnostic(data,'minimum_diagonal');c,d=sector_diagnostic(data,'default')
        self.assertEqual(b['seed'],min((3,6,9,12),key=lambda s:(oracle.action(s).get(s,F(0)),s)))
        self.assertEqual(d['seed'],3)
        self.assertEqual(list(a),list(c))

    def test_target_source_hash_must_match_before_construction(self):
        data,_=fixture()
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory);source=p/'source.json';target=p/'target.json'
            source.write_text(json.dumps(data));target.write_text(json.dumps({'source_sha256':'wrong','target_rotation':[[1,0],[0,1]]}))
            with self.assertRaises(ValueError): construct(source,target,p/'out')
            self.assertFalse((p/'out').exists())


if __name__=='__main__': unittest.main()
