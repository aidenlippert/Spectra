from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import unittest

from experiments.marginal_symbolic import add,adj,mono,scale,product,encode
from experiments.marginal_transfer_verify import apply_word
from research.molecular_collective_20260913.core import (
    density,square_form,extract,tail_replay,rank_replay,support_replay,retained_polynomial,digest)

ROOT=Path(__file__).resolve().parents[2]
CAMPAIGN=ROOT/'results/molecular_collective_20260913/campaign'


def saved(rank=10):
    directory=CAMPAIGN/'h6';data=json.loads((directory/'fixture.json').read_text())
    return data,json.loads((directory/f'rank_{rank}/tail.json').read_text())


class AlgebraTests(unittest.TestCase):
    def test_spin_orbital_casimir_identity(self):
        for s in (2,3):
            n=add(*(density(p,p) for p in range(s)))
            sz=add(*(mono(((1,2*p+t),(0,2*p+t)),F(1-2*t,2)) for p in range(s) for t in range(2)))
            sp=add(*(mono(((1,2*p),(0,2*p+1))) for p in range(s)));sm=adj(sp)
            spin=add(product(sz,sz),scale(add(product(sp,sm),product(sm,sp)),F(1,2)))
            orbital=add(*(product(density(p,q),density(q,p)) for p in range(s) for q in range(s)))
            self.assertFalse(add(orbital,scale(spin,2),scale(n,-s-2),scale(product(n,n),F(1,2))))
            symmetric=add(*(product(density(p,p),density(p,p)) for p in range(s)),
                          *(scale(product(add(density(p,q),density(q,p)),add(density(p,q),density(q,p))),F(1,2)) for p,q in combinations(range(s),2)))
            anti=add(*(scale(product(adj(add(density(p,q),scale(density(q,p),-1))),add(density(p,q),scale(density(q,p),-1))),F(1,2)) for p,q in combinations(range(s),2)))
            self.assertFalse(add(orbital,scale(symmetric,-1),scale(anti,-1)))

    def test_number_centering_is_exact_on_fixed_sector(self):
        pairs=[(i,j) for i in range(3) for j in range(i,3)]
        basis=[density(i,j) if i==j else add(density(i,j),density(j,i)) for i,j in pairs]
        matrix=[[F((i+1)*(j+1),13)+F(i==j,7) for j in range(6)] for i in range(6)]
        h=add(square_form(matrix,basis),mono(((1,0),(0,0)),F(-2,3)))
        for number in (0,1,3,6):
            data={'modes':6,'particles':number,'hamiltonian':encode(h)};p=extract(data)
            self.assertEqual(p['eta'],0)
            restored=add(p['one'],square_form(p['matrix'],p['basis']))
            delta=add(h,scale(restored,-1))
            for inds in combinations(range(6),number):
                state=sum(1<<i for i in inds);action={}
                for w,c in delta.items():
                    target=apply_word(w,state)
                    if target:
                        t,sign=target;action[t]=action.get(t,F(0))+sign*c
                self.assertFalse(any(action.values()))

    def test_input_budget_and_particle_domain(self):
        for m,n in ((22,6),(7,4),(12,13),(True,0)):
            with self.assertRaises(ValueError):extract({'modes':m,'particles':n,'hamiltonian':[]})


@unittest.skipUnless(CAMPAIGN.exists(),'Persisted campaign required')
class CertificateTests(unittest.TestCase):
    def test_saved_h6_remainder_is_below_budget(self):
        data,cert=saved();r=tail_replay(data,cert)
        self.assertLess(F(r['tail_interval_width_Ha']),F(16,10000));self.assertEqual(r['many_body_states_enumerated'],0)

    def test_wrong_fixture_rejected(self):
        data,cert=saved();data['particles']=5
        with self.assertRaisesRegex(ValueError,'binding'):tail_replay(data,cert)

    def test_zero_tail_padding_rejected(self):
        data,cert=saved();cert['alpha']='0';cert['beta']='0'
        with self.assertRaisesRegex(ValueError,'PSD'):tail_replay(data,cert)

    def test_centered_trace_corruption_rejected(self):
        data,cert=saved();cert['factors'][0]['vector'][0]='1'
        with self.assertRaisesRegex(ValueError,'trace'):tail_replay(data,cert)

    def test_negative_factor_rejected(self):
        data,cert=saved();cert['factors'][0]['weight']='-1'
        with self.assertRaisesRegex(ValueError,'Nonnegative'):tail_replay(data,cert)

    def test_out_of_budget_factor_count_rejected(self):
        data,cert=saved();cert['factors']=cert['factors']*3
        with self.assertRaisesRegex(ValueError,'budget'):tail_replay(data,cert)

    def test_rank_obstruction_and_corruption(self):
        data,_=saved();c=json.loads((CAMPAIGN/'h6/rank_obstruction.json').read_text())
        self.assertEqual(rank_replay(data,c)['minimum_factors'],10)
        c['subspace']=[['0']*10 for _ in c['subspace']]
        with self.assertRaisesRegex(ValueError,'Strict'):rank_replay(data,c)

    def test_support_dual_exact_ceiling(self):
        data,tail=saved();cert=json.loads((CAMPAIGN/'h6/rank_10/support.json').read_text())
        r=support_replay(data,tail,cert)
        self.assertGreaterEqual(F(r['support_family_ceiling_Ha']),F(r['certified_lower_Ha']))
        self.assertLess(F(r['support_primal_dual_gap_Ha']),F(1,1000000))

    def test_support_bad_density_rejected(self):
        data,tail=saved();cert=json.loads((CAMPAIGN/'h6/rank_10/support.json').read_text())
        cert['dual_density'][0][0]='3'
        with self.assertRaisesRegex(ValueError,'trace'):support_replay(data,tail,cert)

    def test_support_wrong_tail_rejected(self):
        data,tail=saved();cert=json.loads((CAMPAIGN/'h6/rank_10/support.json').read_text())
        tail['beta']=str(F(tail['beta'])+1)
        with self.assertRaisesRegex(ValueError,'binding'):support_replay(data,tail,cert)

    def test_one_body_false_psd_rejected(self):
        data,tail=saved();cert=json.loads((CAMPAIGN/'h6/rank_10/support.json').read_text())
        cert['chemical_potential']='100'
        with self.assertRaisesRegex(ValueError,'PSD'):support_replay(data,tail,cert)


if __name__=='__main__':unittest.main()
