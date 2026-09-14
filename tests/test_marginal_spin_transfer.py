import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_determinant_tree import DeterminantOracle
from experiments.marginal_fixed_point_ldl import propose
from experiments.marginal_operator_norm import hopping_norm, verify_norm, replay as norm_replay
from experiments.marginal_spin_constructor import wrap
from experiments.marginal_spin_reduction import SpinZeroOracle, replay, spin_operators
from experiments.marginal_spin_transfer import transfer
from experiments.marginal_symbolic import add, mono, encode, scale


def reference_fixture():
    h = add(*(mono(((1,i),(0,i)), 1 if i < 2 else 2) for i in range(4)))
    base = {'modes': 4, 'particles': 2, 'hamiltonian': encode(h)}
    c = dict(base, retained_states=[3], complement_lower='5/2', lower='199/100',
             independent_upper={'states':[3], 'amplitudes':[1]},
             response_basis=[{'states':[6], 'amplitudes':[1]}],
             blocks=[{'states':[s], 'factor':propose([[d-F(5,2)]])} for s,d in [(6,3),(9,3),(12,4)]])
    return wrap(base, c)


def connected_fixture():
    previous = reference_fixture()
    base = previous['spin_symmetric_certificate']
    h = SpinZeroOracle(base).h
    delta = add(*(mono(((1,i),(0,j)), F(1,10)) for i,j in [(0,2),(2,0),(1,3),(3,1)]))
    current = {key: previous[key] for key in ('modes','particles')}
    current['hamiltonian'] = encode(add(h, delta))
    reduced = dict(current, retained_states=[3], lower='39/20', complement_lower='21/10',
        independent_upper={'states':[3], 'amplitudes':[1]},
        response_basis=[{'states':[s], 'amplitudes':[1]} for s in (6,9,12)],
        complement_reference={key:base[key] for key in ('hamiltonian','complement_lower','blocks')})
    return wrap(current, reduced)


class SpinTransferTests(unittest.TestCase):
    def test_two_sided_norm_squares_and_false_identities(self):
        for t in (F(1,10), F(-1,50), F(0)):
            delta = add(*(mono(((1,i),(0,j)), t) for i,j in [(0,2),(2,0),(1,3),(3,1)]))
            certificate = hopping_norm(t)
            self.assertEqual(verify_norm(delta, 4, certificate)[0], 2*abs(t))
            bad = copy.deepcopy(certificate)
            bad['bound'] = str(2*abs(t)+1)
            with self.assertRaises(ValueError): verify_norm(delta, 4, bad)
            bad = copy.deepcopy(certificate)
            bad['minus'][0]['weight'] = '-1'
            with self.assertRaises(ValueError): verify_norm(delta, 4, bad)
        c = connected_fixture()
        reduced = c['spin_symmetric_certificate']
        reduced['complement_lower'] = '23/10'
        with self.assertRaises(ValueError): replay(c)
        reduced['complement_reference']['norm_certificate'] = hopping_norm(F(1,10))
        self.assertEqual(F(replay(c)['complement_coverage']['reference_perturbation_norm_bound']), F(1,5))
        root = Path(__file__).resolve().parents[1]
        c = json.loads((root/'results/marginal_h6/hopping_norm_squares/certificate.json').read_text())
        self.assertEqual(F(norm_replay(c)['operator_norm']), F(1,25))
        c['saturating_witness']['amplitudes'] = [1,0,0,0]
        with self.assertRaises(ValueError): norm_replay(c)

    def test_saved_connected_h6_intervals_keep_small_reference_factors(self):
        root = Path(__file__).resolve().parents[1]
        for name, directions in (('connected_spin',28), ('connected_spin_1_100',31),
                                 ('connected_spin_squares_1_50',31)):
            with self.subTest(name=name):
                path = root/'results/marginal_h6'/name
                c = json.loads((path/'certificate.json').read_text())
                r = replay(c)
                saved = json.loads((path/'receipt.json').read_text())
                self.assertEqual(r['width'], saved['width'])
                self.assertLess(r['width_float'], 3e-10)
                self.assertEqual(r['response_dimension'], directions)
                self.assertEqual(saved['current_q_component_dimensions'], [368])
                self.assertEqual(saved['largest_factor_dimension'], 200)
                self.assertEqual(r['reference_action_states'], 368)
                self.assertEqual(r['unique_determinant_sources'], 400)

    def test_connected_reference_gap_and_separate_original_transfer(self):
        c = connected_fixture()
        r = replay(c)
        self.assertEqual(F(r['width']), F(1,20))
        self.assertEqual(F(r['complement_coverage']['reference_perturbation_norm_bound']), F(2,5))
        self.assertEqual(r['reference_action_states'], 3)
        self.assertEqual(r['unique_determinant_sources'], 4)
        _,_,z = spin_operators(4)
        c['hamiltonian'] = encode(add(DeterminantOracle(c).h, scale(z,F(1,100))))
        r = replay(c)
        self.assertEqual(F(r['spin_symmetry_error_bound']), F(1,50))
        self.assertEqual(F(r['complement_coverage']['reference_perturbation_norm_bound']), F(2,5))
        self.assertEqual(F(r['lower']), F(193,100))

    def test_false_threshold_reference_and_actual_response_are_rejected(self):
        c = connected_fixture()
        for update in ({'complement_lower':'210001/100000'}, {'lower':'199/100'},
                       {'complement_reference':None}, {'retained_states':[3,3]}):
            bad = copy.deepcopy(c)
            bad['spin_symmetric_certificate'].update(update)
            with self.assertRaises(ValueError): replay(bad)
        for update in ({'blocks':[]}, {'complement_lower':'4'},
                       {'hamiltonian':encode(mono(((1,0),(0,0))))}):
            bad = copy.deepcopy(c)
            bad['spin_symmetric_certificate']['complement_reference'].update(update)
            with self.assertRaises(ValueError): replay(bad)
        c['spin_symmetric_certificate']['complement_lower'] = '2'
        self.assertEqual(F(replay(c)['width']), F(1,20))

    def test_transfer_constructs_upper_and_response_without_new_factors(self):
        c = reference_fixture()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'source.json'
            source.write_text(json.dumps(c))
            with patch('experiments.marginal_fixed_point_ldl.propose', side_effect=AssertionError('New factorization')), \
                 contextlib.redirect_stdout(io.StringIO()):
                transfer(source, root/'out', F(1,10), upper_steps=3)
            r = json.loads((root/'out/receipt.json').read_text())
            certificate = json.loads((root/'out/certificate.json').read_text())
            self.assertLess(r['width_float'], 1e-7)
            self.assertEqual(r['current_q_component_dimensions'], [3])
            self.assertEqual(r['largest_factor_dimension'], 1)
            self.assertEqual(r['construction_gap_reference_action_states'], 3)
            self.assertEqual(r['source_verification_action_states']['source_spin_symmetric'], 4)
            self.assertEqual(replay(certificate)['width'], r['width'])
            with self.assertRaises(ValueError): transfer(source, root/'out')
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
                transfer(source, root/'coefficient_reject', F(1,5), upper_steps=3)
            with contextlib.redirect_stdout(io.StringIO()):
                transfer(source, root/'norm_squares', F(1,5), upper_steps=3, norm_squares=True)
            sharp = json.loads((root/'norm_squares/certificate.json').read_text())
            self.assertLess(replay(sharp)['width_float'], 1e-7)


if __name__ == '__main__':
    unittest.main()
