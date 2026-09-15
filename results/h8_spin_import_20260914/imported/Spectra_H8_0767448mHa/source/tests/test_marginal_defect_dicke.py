import copy
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import random
import unittest

from experiments.marginal_defect_dicke import DefectDicke, combine, prepare, replay, enlarged_replay
from experiments.marginal_general_schur import apply_columns, gram, fixture, prepare as explicit_prepare
from experiments.marginal_collective import hopping_polynomial
from experiments.marginal_symbolic import add, mono

ROOT = Path(__file__).resolve().parents[1]


def expand(model, vector):
    """Test oracle only: enumerate spectators and transform pair-order gauge."""
    result = {}
    m = model.pairs
    for (occupation_tuple, k), coefficient in vector.items():
        fixed = dict(occupation_tuple)
        spectators = [i for i in range(m) if i not in fixed]
        for rights in combinations(spectators, k):
            occupations = dict(fixed)
            occupations.update({i: 2 if i in rights else 1 for i in spectators})
            chosen = [mode for i in range(m) for bit, mode in ((1, i), (2, i + m)) if occupations[i] & bit]
            parity = sum(a > b for i, a in enumerate(chosen) for b in chosen[i + 1:])
            state = sum(1 << mode for mode in chosen)
            result[state] = result.get(state, F(0)) + coefficient * (-1) ** parity
    return {state: value for state, value in result.items() if value}


class DefectDickeTests(unittest.TestCase):
    def test_kernel_coordinates_remain_valid_after_growth_and_reordering(self):
        model = DefectDicke(4)
        first = model.atom(rights=2)
        cached = model._coordinate_state(first)
        atoms = [first] + [model.atom(((i, o),), k)
                           for i in range(4) for o in range(4) for k in range(4)]
        expected = gram([expand(model, a) for a in atoms], [expand(model, a) for a in atoms])
        self.assertEqual(model.gram(atoms, atoms), expected)
        self.assertEqual(model._coordinate_inner(cached, cached), 6)
        for a, overlap in zip(atoms, expected[0]):
            self.assertEqual(model._coordinate_inner(cached, model._coordinate_state(a)), overlap)
        reversed_model = DefectDicke(4)
        for a in reversed(atoms):
            reversed_model._coordinate_state(a)
        self.assertEqual(reversed_model.gram(atoms, atoms), expected)

    def test_redundant_kernel_atoms_do_not_grow_physical_rank(self):
        model = DefectDicke(4)
        left = model.atom(((0, 1),), 2)
        right = model.atom(((0, 2),), 1)
        model.inner(left, right)
        rank = sum(len(s['diagonal']) for s in model.coordinate_sectors.values())
        reference = model.atom(rights=2)
        self.assertEqual(model.inner(reference, reference), 6)
        self.assertEqual(sum(len(s['diagonal']) for s in model.coordinate_sectors.values()), rank)
        null = combine((1, reference), (-1, left), (-1, right))
        self.assertEqual(model.inner(null, null), 0)

    def test_overlaps_of_distinct_defect_sets_and_null_relations(self):
        model = DefectDicke(3)
        atoms = [model.atom(rights=k) for k in range(4)]
        atoms += [model.atom(((i, o),), k) for i in range(3) for o in range(4) for k in range(3)]
        explicit = [expand(model, a) for a in atoms]
        self.assertEqual(model.gram(atoms, atoms), gram(explicit, explicit))
        null = combine((1, model.atom(rights=1)), (-1, model.atom(((0, 1),), 1)), (-1, model.atom(((0, 2),), 0)))
        self.assertTrue(null)  # redundant coordinates need not cancel as a dict
        self.assertEqual(model.inner(null, null), 0)
        self.assertFalse(expand(model, null))

    def test_all_single_car_actions_including_empty_and_double_defects(self):
        model = DefectDicke(3)
        atoms = model.reference() + [model.atom(((i, o),), k) for i in range(3) for o in range(4) for k in range(3)]
        for atom in atoms:
            for creation, mode in product((0, 1), range(6)):
                word = ((creation, mode),)
                self.assertEqual(expand(model, model.apply_word(word, atom)), apply_columns(mono(word), [expand(model, atom)])[0])

    def test_cross_pair_car_products_and_exact_anticommutators(self):
        model = DefectDicke(4)
        rng = random.Random(314159)
        for degree in (2, 4, 8):
            for _ in range(35):
                word = tuple((rng.randrange(2), rng.randrange(8)) for _ in range(degree))
                atom = model.atom(((rng.randrange(4), rng.randrange(4)),), rng.randrange(4))
                self.assertEqual(expand(model, model.apply_word(word, atom)), apply_columns(mono(word), [expand(model, atom)])[0])
        for p, q in product(range(8), repeat=2):
            for v in model.reference():
                state = combine((1, model.apply_word(((0, p), (1, q)), v)), (1, model.apply_word(((1, q), (0, p)), v)), (-int(p == q), v))
                self.assertEqual(model.inner(state, state), 0)

    def test_fast_reference_action_on_every_local_occupation_pattern(self):
        model = DefectDicke(4)
        h0 = hopping_polynomial(8, F(1, 5))
        for a, b in product(range(4), repeat=2):
            for k in range(3):
                v = model.atom(((0, a), (2, b)), k)
                expected = apply_columns(h0, [expand(model, v)])[0]
                self.assertEqual(expand(model, model.reference_action(v)), expected)
                other = model.apply_polynomial(h0, v)
                self.assertEqual(model.inner(combine((1, other), (-1, model.reference_action(v))), combine((1, other), (-1, model.reference_action(v)))), 0)

    def test_coupling_matches_explicit_signed_fock_calculation(self):
        model = DefectDicke(5)
        for mixed in (False, True):
            h = fixture(F(1, 100), mixed)
            got = model.coupling(h)
            expected = explicit_prepare(h, 'hermitian_pairs')
            for key in ('metric', 'projected_h', 'perturbation_norm_bound'):
                self.assertEqual(got[key], expected[key])
            self.assertEqual([expand(model, w) for w in got['coupling']], expected['coupling'])
            self.assertEqual(model.gram(got['coupling'], got['coupling']), expected['leakage'])

    def test_reference_action_and_overlap_at_one_hundred_pairs(self):
        from math import comb
        model = DefectDicke(100)
        v = model.atom(rights=50)
        self.assertEqual(model.inner(v, v), comb(100, 50))
        hv = model.reference_action(v)
        self.assertEqual(len(hv), 3)
        self.assertEqual(model.inner(v, hv), 2450 * comb(100, 50))
        word = ((1, 0), (0, 101))
        w = model.apply_word(word, v)
        self.assertLessEqual(len(w), 4)
        self.assertGreater(model.inner(w, w), 0)

    def test_input_and_degree_budgets(self):
        with self.assertRaises(ValueError): DefectDicke(1)
        model = DefectDicke(3)
        with self.assertRaises(ValueError): model.atom(((0, 1), (0, 2)), 1)
        with self.assertRaises(ValueError): model.apply_word(((1, 6),), model.atom(rights=1))
        with self.assertRaises(ValueError): model.coupling(mono(((1, 0),)))
        with self.assertRaises(ValueError): model.recurrence(model.reference(), max_degree=1)
        self.assertEqual(model.recurrence([{}])['annihilator'], [F(1)])

    def test_recurrence_moments_and_saved_certificates_match_explicit(self):
        for name, mixed in (('cycle', False), ('mixed', True)):
            h = fixture(F(1, 1000), mixed)
            got = prepare(h)
            expected = explicit_prepare(h, 'hermitian_pairs', True)
            for key in ('annihilator', 'moments', 'leakage'):
                self.assertEqual(got[key], expected[key])
            certificate = json.loads((ROOT / 'results/marginal_general_schur' / (name + '_1_1000_paired_resolvent') / 'certificate.json').read_text())
            receipt = replay(certificate)
            self.assertLess(receipt['width_float'], 1e-6)
            self.assertEqual(got['zero_residual_norm'], '0')
        changed = copy.deepcopy(certificate)
        changed['lower'] = '4'
        with self.assertRaises(ValueError): replay(changed)
        changed = copy.deepcopy(certificate)
        changed['independent_upper']['amplitudes'] = [1]
        with self.assertRaises(ValueError): replay(changed)

    def test_closure_uses_physical_rank_and_enlarged_certificate_replays(self):
        model = DefectDicke(3)
        v = model.atom(rights=1)
        alternative = combine((1, model.atom(((0, 1),), 1)), (1, model.atom(((0, 2),), 0)))
        basis = model.closure([v, alternative])
        self.assertEqual(len(basis), 4)
        self.assertEqual(len(model.closure([v])), 4)
        with self.assertRaises(ValueError): model.closure([v], max_dimension=3)
        certificate = json.loads((ROOT / 'results/marginal_enlarged_schur/cycle_1_100/certificate.json').read_text())
        receipt = enlarged_replay(certificate)
        self.assertEqual(receipt['retained_dimension'], 14)
        self.assertLess(receipt['width_float'], 3e-5)
        changed = copy.deepcopy(certificate)
        changed['lower'] = '4'
        with self.assertRaises(ValueError): enlarged_replay(changed)
        changed = copy.deepcopy(certificate)
        changed['rounds'] = 2
        with self.assertRaises(ValueError): enlarged_replay(changed)


if __name__ == '__main__':
    unittest.main()
