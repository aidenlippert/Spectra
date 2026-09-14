import json
import unittest
from fractions import Fraction as F
from pathlib import Path

from experiments.marginal_hopping_filter import verify_edge_identity, merge_shift, replay
from experiments.marginal_transfer_verify import apply_word


ROOT = Path(__file__).resolve().parents[1]
VALENCE = {90:-1,102:2,105:-1,150:-1,153:2,165:-1}


def hopping(cut,coefficient=-1):
    terms = []
    for spin in (0,1):
        left,right = 2*(cut-1)+spin,2*cut+spin
        terms += [(((1,left),(0,right)),coefficient),
                  (((1,right),(0,left)),coefficient)]
    return terms


def action(terms,state):
    out = {}
    for s,a in state.items():
        for word,b in terms:
            image = apply_word(word,s)
            if image:
                t,phase = image
                out[t] = out.get(t,0)+a*b*phase
    return {s:a for s,a in out.items() if a}


def dot(a,b):
    return sum(x*b.get(s,0) for s,x in a.items())


class HoppingFilterTests(unittest.TestCase):
    def test_edge_identity_exact_three_sizes(self):
        for sites in (2,4,8):
            for edge in (0,sites-1):
                for spin in (0,1):
                    self.assertEqual(verify_edge_identity(sites,edge,spin)['residual_terms'],0)
        for args in ((12,0,0),(4,1,0),(4,0,2)):
            with self.assertRaises(ValueError):verify_edge_identity(*args)

    def test_merge_shift_is_exact(self):
        self.assertEqual(merge_shift(4,'1/5'),F(-3,13))
        self.assertEqual(merge_shift(0,'1/5'),F(-5,13))
        with self.assertRaises(ValueError):merge_shift(4,.2)

    def test_two_and_three_correlated_blocks_by_direct_car(self):
        # Explicit sparse finite states independently test the induction after
        # a cluster has become internally entangled. No global sector matrix.
        for q in (2,3):
            with self.subTest(blocks=q):
                sites = 4*q
                state = {0:1}
                for block in range(q):
                    state = {s+(t<<(8*block)):a*b for s,a in state.items()
                             for t,b in VALENCE.items()}
                for cut in range(4,sites,4):
                    moved = action(hopping(cut),state)
                    state = {s:5*state.get(s,0)-moved.get(s,0) for s in state.keys()|moved.keys()}
                    state = {s:a for s,a in state.items() if a}
                    self.assertLessEqual(len(state),65536)
                norm = dot(state,state)
                self.assertEqual(norm,12**q*26**(q-1))
                self.assertTrue(all(s.bit_count()==sites for s in state))
                self.assertTrue(all(sum((s>>(2*i))&1 for i in range(sites))==sites//2 for s in state))
                # Full Hubbard energy: onsite cost plus ALL internal hoppings.
                onsite = sum(a*a*4*sum(((s>>(2*i))&3)==3 for i in range(sites))
                             for s,a in state.items())
                kinetic = sum(dot(state,action(hopping(cut),state)) for cut in range(1,sites))
                self.assertEqual(F(onsite+kinetic,norm),-(q-1)*F(3,13))
                # Remote endpoint observables are inherited exactly.
                for edge in (0,sites-1):
                    for spin in (0,1):
                        occupation = sum(a*a*((s>>(2*edge+spin))&1) for s,a in state.items())
                        self.assertEqual(F(occupation,norm),F(1,2))
                    doublons = sum(a*a*int(((s>>(2*edge))&3)==3) for s,a in state.items())
                    self.assertEqual(doublons,0)
                for cut in (1,sites-1):
                    self.assertEqual(dot(state,action(hopping(cut,1),state)),0)

    def test_direct_merge_validation_uses_a_noneigenstate(self):
        # The boundary state is NOT an eigenstate: its energy is zero, while
        # the full H action has nonzero charge fluctuations.
        hv = action(sum((hopping(cut) for cut in range(1,4)),[]),VALENCE)
        self.assertGreater(dot(hv,hv),0)
        self.assertEqual(dot(VALENCE,hv),0)

    def test_actual_h8_correlated_chain_and_refusals(self):
        c = json.loads((ROOT/'results/marginal_graded_hubbard8/tiled_eight_upper/certificate.json').read_text())
        result = replay(c['upper'],c['hamiltonian'],'23/100',1000000)
        e8 = F(result['boundary_data']['energy8'])
        self.assertLess(F(result['upper_per_site']),e8/8-F(28,1000))
        self.assertEqual(result['filter_norm_exponent'],124999)
        self.assertEqual(result['boundary_data']['edge_spin_occupations'],[['1/2','1/2']]*2)
        self.assertLessEqual(result['boundary_data']['orbit_support'],4096)
        with self.assertRaises(ValueError):replay(c['upper'],c['hamiltonian'],.23,16)
        with self.assertRaises(ValueError):replay(c['upper'],c['hamiltonian'],'23/100',15)
        with self.assertRaises(ValueError):replay(c['upper'],c['hamiltonian'],'23/100',18)


if __name__ == '__main__':
    unittest.main()
