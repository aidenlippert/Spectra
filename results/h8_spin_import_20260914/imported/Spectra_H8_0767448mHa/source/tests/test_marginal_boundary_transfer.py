import unittest
from fractions import Fraction as F

from experiments.marginal_boundary_transfer import compile_state, contract
from experiments.marginal_local_hubbard_block import _terms
from experiments.marginal_transfer_verify import apply_word


VALENCE = {0x99: 1, 0x96: -1, 0x69: -1, 0x66: 1}


def act(terms, state):
    out = {}
    for s, amplitude in state.items():
        for word, coefficient in terms:
            image = apply_word(word, s)
            if image:
                target, sign = image
                out[target] = out.get(target, 0) + amplitude * coefficient * sign
    return {s: a for s, a in out.items() if a}


def blocks_state(blocks,seed=VALENCE):
    state = {0: 1}
    for block in range(blocks):
        state = {s | (t << (8 * block)): a * b
                 for s, a in state.items() for t, b in seed.items()}
    return state


def filtered_state(state, blocks, a, b):
    for cut in range(1, blocks):
        terms = []
        for spin in (0, 1):
            left, right = 2 * (4 * cut - 1) + spin, 2 * (4 * cut) + spin
            terms.extend([(((1, left), (0, right)), -1),
                          (((1, right), (0, left)), -1)])
        first = act(terms, state)
        second = act(terms, first)
        keys = state.keys() | first.keys() | second.keys()
        state = {s: state.get(s, 0) - a * first.get(s, 0) + b * second.get(s, 0)
                 for s in keys}
        state = {s: v for s, v in state.items() if v}
    return state


def reference(blocks, a, b, seed=VALENCE):
    sites = 4 * blocks
    state = filtered_state(blocks_state(blocks,seed), blocks, F(a), F(b))
    if len(state)>100000:raise ValueError('Direct test support cap exceeded')
    norm = sum(x * x for x in state.values())
    onsite = [4] * sites
    hopping = [1] * (sites - 1)
    h = _terms(sites, 4, 1, onsite, hopping)
    moved = act(h, state)
    numerator = sum(x * moved.get(s, 0) for s, x in state.items())
    return norm, numerator


class BoundaryTransferTests(unittest.TestCase):
    def test_direct_car_reference_for_one_two_three_blocks(self):
        for a, b in ((0, 0), (F(1, 5), 0), (F(1, 5), F(1, 7)), (0, 1)):
            for blocks in (1, 2, 3):
                with self.subTest(a=a, b=b, blocks=blocks):
                    state = blocks_state(blocks)
                    self.assertLessEqual(len(state), 100000)
                    compiled = compile_state(blocks_state(1), 4)
                    actual = contract(compiled, a, b, blocks)
                    norm, numerator = reference(blocks, a, b)
                    self.assertEqual(F(actual['norm']), F(norm))
                    self.assertEqual(F(actual['energy_numerator']), F(numerator))
                    self.assertEqual(actual['energy'], str(F(numerator, norm)))

    def test_charged_noneigenstate_by_direct_car(self):
        moved=act(_terms(4,4,1),VALENCE)
        seed={s:8*VALENCE.get(s,0)-moved.get(s,0) for s in VALENCE.keys()|moved.keys()}
        seed={s:a for s,a in seed.items() if a}
        self.assertTrue(any((s&3)==3 for s in seed))
        compiled=compile_state(seed,4)
        for q in (1,2,3):
            actual=contract(compiled,'1/5','1/7',q)
            norm,numerator=reference(q,'1/5','1/7',seed)
            self.assertEqual(F(actual['norm']),norm)
            self.assertEqual(F(actual['energy_numerator']),numerator)

    def test_refusal_gates(self):
        compiled = compile_state(blocks_state(1), 4)
        with self.assertRaises(ValueError): compile_state({0: 1}, 4)
        with self.assertRaises(ValueError): compile_state(blocks_state(1), 6)
        with self.assertRaises(ValueError): contract(compiled, 0.2, 0, 1)
        with self.assertRaises(ValueError): contract(compiled, 0, 0.1, 1)
        with self.assertRaises(ValueError): contract(compiled, 0, 0, 0)
        with self.assertRaises(ValueError): contract(compiled, 0, 0, 65)


if __name__ == '__main__':
    unittest.main()
