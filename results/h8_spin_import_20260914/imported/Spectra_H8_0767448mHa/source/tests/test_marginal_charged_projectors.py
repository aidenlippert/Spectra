import unittest
from fractions import Fraction as F

from experiments.marginal_charged_projectors import (
    charged_vectors, charged_overlap_grams, charged_projector_bound, _particlehole,
    joint_projector_bound,
)
from experiments.marginal_projector_extendibility import _psd
from experiments.marginal_projector_extendibility import replay


SOURCE = {62: 1, 3008: 1}


def spinflip(s):
    target=0;sign=1
    for site in range(6):
        bits=(s>>(2*site))&3
        target|=((bits&1)*2+(bits>>1))<<(2*site)
        if bits==3:sign=-sign
    return target,sign


def particlehole(s):
    target=4095;amplitude=1
    for mode in range(11,-1,-1):
        if s>>mode&1:
            amplitude*=(-1)**(mode//2)
            amplitude*=(-1)**((target&((1<<mode)-1)).bit_count())
            target^=1<<mode
    return target,amplitude


def transform(vector, fn):
    return {fn(s)[0]: fn(s)[1] * a for s, a in vector.items()}


def independent_vectors():
    first = SOURCE
    second = transform(first, spinflip)
    third = transform(first, particlehole)
    fourth = transform(second, particlehole)
    return [first, second, third, fourth]


def independent_grams(windows):
    assert windows==2
    vectors = independent_vectors()
    groups = {}
    for offset in range(2):
        for vector in vectors:
            for environment in range(4):
                column = ({s+4096*environment:a for s,a in vector.items()} if offset==0
                          else {4*s+environment:a for s,a in vector.items()})
                # Full spin-sector key, computed directly from seven-site bits.
                sectors = set()
                for s in column:
                    up = sum((s >> (2 * i)) & 1 for i in range(windows + 5))
                    sectors.add((up, s.bit_count() - up))
                assert len(sectors) == 1
                groups.setdefault(next(iter(sectors)), []).append(column)
    return {key: [[sum(a * right.get(s, 0) for s, a in left.items())
                    for right in columns] for left in columns]
            for key, columns in sorted(groups.items())}


class ChargedProjectorTests(unittest.TestCase):
    def joint_certificate(self):
        return {'kind':'hubbard_projector_extension_v4','chain_sites':10,
                'target':{'U':'0','t':'0','V':'0'},
                'local_window':{'kind':'local_hubbard_block_v1','sites':6,'U':'0','t':'0','V':'0'},
                'vector':{0x999:1,0x666:1},'windows':2,'projector_sum_ceiling':'2',
                'penalty':'1/3','penalized_lower':'0',
                'joint':{'vector':SOURCE,'windows':2,'ratio':'1/2',
                         'projector_sum_ceiling':'2','penalty':'1/5'}}

    def test_joint_energy_full_fock_and_translation_accounting(self):
        # K=0 gives a known all-Fock PSD operator: a sum of orthogonal
        # positive projectors. The vacuum fixes its exact minimum to zero.
        result=replay(self.joint_certificate())
        self.assertTrue(result['accepted'])
        self.assertEqual(result['local_sum_dimensions'],4096)
        self.assertEqual(len(result['local_sectors']),94)
        self.assertEqual(result['local_maximum_psd_dimension'],200)
        self.assertEqual(F(result['half_local_penalty']),F(8,15))
        self.assertEqual(F(result['charged_local_penalty']),F(1,10))
        self.assertEqual(F(result['periodic_lower_density']),-F(8,75))
        self.assertEqual(F(result['open_lower_energy']),-F(16,15))

    def test_joint_energy_refuses_invalid_proofs(self):
        for update in ({'penalty':'-1'}, {'ratio':'0'}, {'ratio':'-1'},
                       {'windows':5}, {'vector':{62:1}}, {'projector_sum_ceiling':'1'}):
            c=self.joint_certificate();c['joint'].update(update)
            with self.subTest(update=update),self.assertRaises(ValueError):replay(c)
        for update in ({'penalty':'-1'}, {'penalized_lower':'1/100'},
                       {'vector':{0x999:1}}, {'vector':{63:1}},
                       {'joint':None}, {'kind':'hubbard_projector_extension_v3'}):
            c=self.joint_certificate();c.update(update)
            with self.subTest(update=update),self.assertRaises(ValueError):replay(c)
        c=self.joint_certificate();c['chain_sites']=8;c['joint']['windows']=4
        c['joint']['projector_sum_ceiling']='4'
        with self.assertRaisesRegex(ValueError,'support'):replay(c)
        c=self.joint_certificate();c['target']['U']='4'
        with self.assertRaisesRegex(ValueError,'target'):replay(c)

    def test_particle_hole_matches_direct_car_on_all_fock_states(self):
        for s in range(4096):
            self.assertEqual(_particlehole(s),particlehole(s))

    def test_transforms_and_odd_window_gram_match(self):
        expected_vectors = independent_vectors()
        vectors, norm = charged_vectors(SOURCE)
        self.assertEqual(norm, 2)
        self.assertEqual(vectors, expected_vectors)
        for actual, expected in zip(vectors, expected_vectors):
            self.assertEqual(sum(x * x for x in actual.values()), 2)
            self.assertEqual(actual, expected)
        for i,left in enumerate(vectors):
            for right in vectors[i+1:]:self.assertEqual(sum(a*right.get(s,0) for s,a in left.items()),0)
        expected = independent_grams(2)
        actual_norm, actual = charged_overlap_grams(SOURCE, 2)
        self.assertEqual(actual_norm, 2)
        self.assertEqual(actual, expected)

    def test_full_four_source_window_dimensions(self):
        result = charged_projector_bound(SOURCE, 4, 4)
        self.assertTrue(result['accepted'])
        self.assertEqual(result['gram_dimension'], 1024)
        self.assertEqual(result['maximum_psd_dimension'], 96)
        self.assertEqual(result['support_sites'], 9)

    def test_small_ceiling_acceptance_and_refusals(self):
        result = charged_projector_bound(SOURCE, 2, 2)
        self.assertEqual(result['average_fidelity_ceiling'], '1')
        with self.assertRaises(ValueError): charged_overlap_grams({1: 1}, 2)
        with self.assertRaises(ValueError): charged_overlap_grams({62: 1, 63: 1}, 2)
        with self.assertRaises(ValueError): charged_overlap_grams(SOURCE, 1)
        with self.assertRaises(ValueError): charged_overlap_grams(SOURCE, 5)
        with self.assertRaises(ValueError): charged_projector_bound(SOURCE, 2, 0)
        with self.assertRaises(ValueError): charged_projector_bound(SOURCE, 2, 3)
        with self.assertRaises(ValueError): charged_overlap_grams({62: F(1)}, 2)
        with self.assertRaisesRegex(ValueError,'reflection'):charged_vectors({62:1})

    def test_joint_bound_matches_direct_outer_product(self):
        half = {0x999: 1, 0x666: 2}; families = [half] + independent_vectors(); groups = {}
        norms=[sum(a*a for a in v.values()) for v in families]
        for offset in range(2):
            for index, vector in enumerate(families):
                for environment in range(4):
                    column = ({s+4096*environment:a for s,a in vector.items()} if offset == 0
                              else {4*s+environment:a for s,a in vector.items()})
                    s = next(iter(column)); up = sum((s>>(2*i))&1 for i in range(7))
                    groups.setdefault((up,s.bit_count()-up), []).append((column,index))
        def accepts(fn):
            try:fn()
            except ValueError:return False
            return True
        outcomes=[]
        for ceiling in (F(1), F(3,2), F(2)):
            def direct():
                for cols in groups.values():
                    support=sorted(set().union(*(set(c) for c,_ in cols)))
                    _psd([[ceiling*(x==y)-sum((F(1) if i==0 else F(1,2))*c.get(x,0)*c.get(y,0)/norms[i] for c,i in cols)
                           for y in support] for x in support])
            expected=accepts(direct);outcomes.append(expected)
            self.assertEqual(accepts(lambda:joint_projector_bound(half,SOURCE,2,F(1,2),ceiling)),expected)
            self.assertEqual(accepts(lambda:joint_projector_bound({s:7*a for s,a in half.items()},
                                      {s:13*a for s,a in SOURCE.items()},2,F(1,2),ceiling)),expected)
        self.assertIn(False,outcomes)
        self.assertIn(True,outcomes)


if __name__ == '__main__':
    unittest.main()
