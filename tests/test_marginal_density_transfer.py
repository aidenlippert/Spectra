import unittest
from fractions import Fraction as F

from experiments.marginal_boundary_transfer import compile_state, contract, replay_target
from experiments.marginal_local_hubbard_block import _actions, _terms, replay as replay_local
from experiments.marginal_transfer_verify import apply_word

V4 = {0x99: 1, 0x96: -1, 0x69: -1, 0x66: 1}


def act(terms, state):
    out = {}
    for s, x in state.items():
        for word, a in terms:
            image = apply_word(word, s)
            if image:
                t, sign = image
                out[t] = out.get(t, 0) + x * a * sign
    return {s: x for s, x in out.items() if x}


def blocks(q,seed=V4):
    state = {0: 1}
    for k in range(q):
        state = {s | (x << (8 * k)): a * b
                 for s, a in state.items() for x, b in seed.items()}
    return state


def reference(q, a, b, U, t, V, seed=V4):
    state = blocks(q,seed)
    for cut in range(1, q):
        terms = []
        for spin in (0, 1):
            x, y = 2 * (4 * cut - 1) + spin, 2 * (4 * cut) + spin
            terms += [(((1, x), (0, y)), -1), (((1, y), (0, x)), -1)]
        h1 = act(terms, state)
        h2 = act(terms, h1)
        state = {s: state.get(s, 0) - a * h1.get(s, 0) + b * h2.get(s, 0)
                 for s in state.keys() | h1.keys() | h2.keys()}
        state = {s: x for s, x in state.items() if x}
    sites = 4 * q
    moved = act(_terms(sites, U, t), state)
    numerator = sum(x * moved.get(s, 0) for s, x in state.items())
    for s, x in state.items():
        density = sum((((s >> (2 * i)) & 3).bit_count() - 1) *
                      (((s >> (2 * (i + 1))) & 3).bit_count() - 1)
                      for i in range(sites - 1))
        numerator += x * x * V * density
    return sum(x * x for x in state.values()), numerator


class DensityTransferTests(unittest.TestCase):
    def test_direct_car_density_energies(self):
        moved=act(_terms(4,4,1),V4)
        seed={s:8*V4.get(s,0)-moved.get(s,0) for s in V4.keys()|moved.keys()}
        seed={s:a for s,a in seed.items() if a}
        self.assertTrue(any((s&3)==3 for s in seed))
        for U, t, V in ((4, 1, F(1, 2)), (4, 1, F(-1, 2)),
                        (3, F(2, 3), F(1, 2))):
            for q in (1, 2, 3):
                with self.subTest(U=U, t=t, V=V, blocks=q):
                    actual = contract(compile_state(seed, 4, U, t, V), F(1, 5), F(1, 7), q)
                    norm, numerator = reference(q, F(1, 5), F(1, 7), F(U), F(t), V,seed)
                    self.assertEqual(actual['target'],{'U':str(F(U)),'t':str(F(t)),'V':str(V)})
                    self.assertEqual(F(actual['norm']), norm)
                    self.assertEqual(F(actual['energy_numerator']), numerator)
                    self.assertEqual(actual['energy'], str(F(numerator, norm)))

    def test_local_density_action_diagonal_and_profile_refusal(self):
        for V in (F(1, 2), F(-1, 2)):
            actions = _actions(2, 4, 1, V=V)
            for state, image in actions.items():
                n0 = ((state & 3).bit_count() - 1)
                n1 = (((state >> 2) & 3).bit_count() - 1)
                expected = 4 * (int((state & 3) == 3)+int(((state>>2)&3)==3)) + V * n0 * n1
                expected -= F(4, 2) * (state.bit_count() - 2)
                self.assertEqual(image.get(state, 0), expected)
        with self.assertRaises(ValueError):
            _actions(2, 4, 1, V=1, density=[1, 2])
        with self.assertRaises(ValueError):
            _actions(4,4,1,V=1,density=[0,1,2])

    def test_density_only_all_fock_lower_cannot_ignore_interaction(self):
        for V in (-1,1):
            certificate={'kind':'local_hubbard_block_v1','sites':2,'U':0,'t':0,'V':V,'lower':-1}
            result=replay_local(certificate)
            self.assertEqual(result['sum_dimensions'],16)
            self.assertEqual(result['V'],str(V))
            with self.assertRaises(ValueError):replay_local(dict(certificate,lower='-1/2'))
        certificate.update(V=1,upper_vector={'3':1})
        self.assertEqual(replay_local(certificate)['upper_quotient'],'-1')

    def test_target_window_mismatch_rejected_before_source_work(self):
        window={'kind':'local_hubbard_block_v1','sites':2,'U':2,'t':1,'V':'1/2','lower':-2}
        for changes in ({'U':3},{'t':'2/3'},{'V':'-1/2'}):
            with self.assertRaises(ValueError):
                replay_target({}, {}, 0,0,1,dict(window,**changes),U=4,t=1,V='1/2')


if __name__ == '__main__':
    unittest.main()
