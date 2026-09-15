from fractions import Fraction as F
import copy
import unittest
import json
from pathlib import Path

from experiments.marginal_orbit_certificate import verify, permuted_canonical_word, normalize_scalar_residual
from experiments.marginal_symbolic import add, adj, encode, mono, product, scale, transform


class OrbitCertificateTests(unittest.TestCase):
    def certificate(self,p,perms):
        square=product(adj(p),p)
        h=scale(add(*(transform(square,g) for g in perms)),F(1,len(perms)))
        return {'modes':3,'particles':1,'operator_degree':3,'hamiltonian':encode(h),
                'number_multiplier':[],'b':'0','permutations':perms,
                'orbit_squares':[{'polynomial':encode(p),'weight':'1'}]}

    def test_charged_square_and_fermionic_permutation(self):
        p=add(mono(((0,0),(0,1))),scale(mono(((0,1),(0,2))),F(2)))
        c=self.certificate(p,[[0,1,2],[2,1,0]])
        r=verify(c)
        self.assertEqual(F(r['lower']),0)
        self.assertEqual(F(r['residual_l1']),0)

    def test_single_annihilator_square_is_allowed(self):
        c=self.certificate(mono(((0,0),)),[[1,0,2]])
        self.assertEqual(F(verify(c)['lower']),0)

    def test_asymmetric_direct_squares_are_not_averaged(self):
        c=self.certificate(mono(((0,0),)),[[0,1,2],[1,0,2]])
        p=mono(((0,0),));weight=F(3,7)
        c['hamiltonian']=encode(add(dict(),scale(add(mono(((1,0),(0,0))),mono(((1,1),(0,1)))),F(1,2)),scale(product(adj(p),p),weight)))
        c['direct_squares']=[{'polynomial':encode(p),'weight':str(weight)}]
        self.assertEqual(F(verify(c)['residual_l1']),0)
        bad=copy.deepcopy(c);bad['direct_squares'][0]['weight']='-1'
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(c);bad['direct_squares'][0]['polynomial']=encode(add(mono(()),p))
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(c);bad['direct_squares']={}
        with self.assertRaises(ValueError):verify(bad)

    def test_invalid_charge_weights_and_permutations_rejected(self):
        c=self.certificate(mono(((0,0),)),[[0,1,2]])
        bad=copy.deepcopy(c);bad['orbit_squares'][0]['weight']='-1'
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(c);bad['permutations']=[[False,1,2]]
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(c);bad['permutations']=[[0,0,2]]
        with self.assertRaises(ValueError):verify(bad)
        bad=copy.deepcopy(c);bad['orbit_squares'][0]['polynomial']=encode(add(mono(()),mono(((0,0),))))
        with self.assertRaises(ValueError):verify(bad)

    def test_residual_is_charged_exactly(self):
        c=self.certificate(mono(((0,0),)),[[0,1,2]])
        c['b']='1/10'
        self.assertEqual(F(verify(c)['lower']),0)

    def test_scalar_residual_normalization_preserves_lower_bound(self):
        c=self.certificate(mono(((0,0),)),[[0,1,2]])
        c['b']='100'
        normalized,receipt=normalize_scalar_residual(c)
        self.assertEqual(c['b'],'100')
        self.assertEqual(normalized['b'],'0')
        self.assertEqual(F(receipt['residual_l1']),0)
        self.assertEqual(receipt['lower'],verify(c)['lower'])

    def test_fast_permutation_matches_general_car_normal_ordering(self):
        from itertools import combinations, permutations
        for nc in range(4):
            for na in range(4):
                for left in combinations(range(3),nc):
                    for right in combinations(range(3),na):
                        word=tuple((1,i) for i in left)+tuple((0,i) for i in right)
                        for mapping in permutations(range(3)):
                            image,sign=permuted_canonical_word(word,mapping)
                            self.assertEqual({image:F(sign)},transform(mono(word),mapping))

    def test_saved_compact_certificate_replays(self):
        root=Path(__file__).resolve().parents[1]/'results/marginal_distill'
        certificate=json.loads((root/'compact_certificate.json').read_text())
        saved=json.loads((root/'compact_receipt.json').read_text())
        self.assertEqual(verify(certificate)['lower'],saved['lower'])
        self.assertEqual(len(certificate['orbit_squares']),36)


if __name__=='__main__':unittest.main()
