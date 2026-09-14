import contextlib
import copy
from fractions import Fraction as F
import io
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from experiments.marginal_fixed_point_ldl import propose
from experiments.marginal_spin_temple import inertia, replay, compress
from experiments.marginal_spin_transfer import transfer
from experiments.marginal_symbolic import add, mono, encode


def fixture():
    h = add(*(mono(((1,i),(0,i)), 1 if i < 2 else 2) for i in range(4)))
    base = {'modes':4, 'particles':2, 'hamiltonian':encode(h)}
    reduced = dict(base, retained_states=[3,6,9], complement_lower='7/2',
        excitation_lower='3', response_basis=[],
        independent_upper={'states':[3,12], 'amplitudes':[10,1]},
        blocks=[{'states':[12], 'factor':propose([[F(1,2)]])}])
    return dict(base, kind='spin_temple_interval_v1', spin_symmetric_certificate=reduced)


class SpinTempleTests(unittest.TestCase):
    def test_direct_and_transferred_atom_gap_rebuild_physical_hamiltonian(self):
        c=fixture();reduced=c['spin_symmetric_certificate'];reduced.pop('blocks')
        atoms={'kind':'spin_rational_atom_complement_v1',
               'blocks':[{'states':[12],'metric_weights':[1],'atom_scale':1,'atoms':[]}]}
        reduced['complement_atoms']=atoms
        r=replay(c);self.assertEqual(r['complement_coverage']['blocks'][0]['packing_threshold_lower'],'4')
        self.assertEqual(r['width'],replay(fixture())['width'])
        for change in ({'complement_lower':'5'}, {'blocks':[]}, {'complement_atoms':{'kind':'unknown'}},
                       {'complement_atoms':dict(atoms,blocks=[dict(atoms['blocks'][0],states=[3])])}):
            bad=copy.deepcopy(c);bad['spin_symmetric_certificate'].update(change)
            with self.assertRaises(ValueError): replay(bad)
        transferred=copy.deepcopy(c);red=transferred['spin_symmetric_certificate']
        red.pop('complement_atoms')
        red['complement_reference']={'hamiltonian':red['hamiltonian'],'complement_lower':'7/2','complement_atoms':atoms}
        self.assertEqual(replay(transferred)['width'],r['width'])
        red['complement_lower']='15/4'
        with self.assertRaises(ValueError): replay(transferred)

    def test_exact_inertia_including_zero_pivots_and_nullspaces(self):
        self.assertEqual(inertia([[0,1,0],[1,0,0],[0,0,0]]),
                         {'positive':1, 'negative':1, 'zero':1})
        self.assertEqual(inertia([[0,0],[0,0]]), {'positive':0, 'negative':0, 'zero':2})
        self.assertEqual(inertia([[1,1],[1,1]]), {'positive':1, 'negative':0, 'zero':1})
        rng = random.Random(481)
        for size in range(1,9):
            for _ in range(4):
                diagonal = [rng.choice([-3,-1,0,0,2,4]) for _ in range(size)]
                t = [[(1 if i == j else rng.randint(-3,3) if i < j else 0)
                      for j in range(size)] for i in range(size)]
                matrix = [[sum(t[k][i]*diagonal[k]*t[k][j] for k in range(size))
                           for j in range(size)] for i in range(size)]
                self.assertEqual(inertia(matrix), {'positive':sum(x>0 for x in diagonal),
                    'negative':sum(x<0 for x in diagonal), 'zero':sum(x==0 for x in diagonal)})
        for bad in ([], [[1,2]], [[1,2],[3,1]], [[0]*33 for _ in range(33)]):
            with self.assertRaises(ValueError): inertia(bad)

    def test_temple_with_degenerate_threshold_and_full_witness_variance(self):
        c = fixture()
        r = replay(c)
        mu, variance = F(204,101), F(400,10201)
        self.assertEqual(F(r['spin_symmetric_upper']), mu)
        self.assertEqual(F(r['residual_norm_squared']), variance)
        self.assertEqual(F(r['lower']), mu-variance/(3-mu))
        self.assertLessEqual(F(r['lower']), 2)
        self.assertGreaterEqual(F(r['upper']), 2)
        self.assertEqual(r['schur_lower_inertia'], {'positive':0, 'negative':1, 'zero':2})
        self.assertEqual(r['response_dimension'], 0)
        self.assertEqual(r['witness_support'], 2)  # Includes a state outside P.
        c['spin_symmetric_certificate']['independent_upper']['amplitudes'] = [30,3]
        self.assertEqual(replay(c)['width'], r['width'])
        c['spin_symmetric_certificate']['independent_upper'] = {'states':[3], 'amplitudes':[1]}
        self.assertEqual(F(replay(c)['width']), 0)

    def test_rejects_bad_gap_count_sector_symmetry_and_excited_witness(self):
        for update in ({'excitation_lower':'7/2'}, {'excitation_lower':'31/10'},
                       {'blocks':[]}, {'modes':6}, {'response_basis':None},
                       {'independent_upper':{'states':[6], 'amplitudes':[1]}},
                       {'independent_upper':{'states':[5], 'amplitudes':[1]}},
                       {'hamiltonian':encode(mono(((1,0),(0,0))))}):
            with self.subTest(update=update):
                c = fixture()
                c['spin_symmetric_certificate'].update(update)
                with self.assertRaises(ValueError): replay(c)
        c = fixture()
        mu = F(replay(c)['spin_symmetric_upper'])
        c['spin_symmetric_certificate']['excitation_lower'] = str(mu)
        with self.assertRaises(ValueError): replay(c)
        c['spin_symmetric_certificate']['excitation_lower'] = str(mu+F(1,10**20))
        self.assertLess(F(replay(c)['lower']), -10**10)

    def test_original_hamiltonian_transfer_is_separate(self):
        c = fixture()
        h = add(*(mono(((1,i),(0,i)), 1 if i < 2 else 2) for i in range(4)),
                mono(((1,0),(0,0)), F(1,100)))
        c['hamiltonian'] = encode(h)
        r = replay(c)
        self.assertEqual(F(r['spin_symmetry_error_bound']), F(1,100))
        self.assertEqual(F(r['upper']), F(205,101))
        self.assertEqual(F(r['lower']), F(r['spin_symmetric_lower'])-F(1,100))
        self.assertEqual(F(r['spin_sector_excitation_lower']), 3)

    def test_compression_replays_premises_without_trusting_old_endpoints(self):
        c = fixture()
        c['kind'] = 'spin_reduced_perturbation_interval_v1'
        reduced = c['spin_symmetric_certificate']
        reduced['response_basis'] = [{'states':[12], 'amplitudes':[1]}]
        reduced['lower'] = '999999'  # A proposal source, never an accepted old interval.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'source.json'
            source.write_text(json.dumps(c))
            with contextlib.redirect_stdout(io.StringIO()):
                compress(source, root/'result', F(1,10))
            result = json.loads((root/'result/certificate.json').read_text())
            r = replay(result)
            self.assertEqual(r['response_dimension'], 0)
            self.assertNotIn('lower', result['spin_symmetric_certificate'])
            self.assertLessEqual(F(r['lower']), 2)
            with self.assertRaises(ValueError): compress(source, root/'result')

    def test_saved_h6_prefixes_and_combined_reference_norm_proof(self):
        root = Path(__file__).resolve().parents[1]/'results/marginal_h6'
        for name, count in (('connected_spin_1_100_temple',27), ('connected_spin_squares_1_50_temple',28),
                            ('connected_spin_squares_1_50_temple_fresh',29)):
            with self.subTest(name=name):
                path = root/name
                c = json.loads((path/'certificate.json').read_text())
                saved = json.loads((path/'receipt.json').read_text())
                r = replay(c)
                self.assertEqual(r['width'], saved['width'])
                self.assertEqual(r['response_dimension'], count)
                self.assertEqual(r['schur_lower_inertia'], {'positive':31, 'negative':1, 'zero':0})
                self.assertLess(r['width_float'], 5.5e-11)
                self.assertLess(r['spin_symmetric_width_float'], 1.3e-13)
                self.assertEqual(r['reference_action_states'], 368)
                self.assertEqual(r['unique_determinant_sources'], 400)
        history = json.loads((path/'response_history.json').read_text())
        self.assertEqual(history[0]['response_dimension'], 0)
        self.assertEqual(history[-1]['response_dimension'], 29)
        self.assertEqual(saved['largest_factor_dimension'], 200)
        self.assertEqual(saved['current_q_component_dimensions'], [368])
        self.assertFalse((path/'response_recipe').exists())
        # The larger perturbation needs the sharper norm, not a trusted receipt.
        del c['spin_symmetric_certificate']['complement_reference']['norm_certificate']
        with self.assertRaises(ValueError): replay(c)

    def test_fresh_transfer_does_not_call_full_energy_response_or_new_factors(self):
        from test_marginal_spin_transfer import reference_fixture
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root/'source.json'
            source.write_text(json.dumps(reference_fixture()))
            with patch('experiments.marginal_spin_transfer.refine_response', side_effect=AssertionError('Full response discovery')), \
                 patch('experiments.marginal_fixed_point_ldl.propose', side_effect=AssertionError('New factors')), \
                 contextlib.redirect_stdout(io.StringIO()):
                transfer(source, root/'out', F(1,5), upper_steps=3, norm_squares=True, temple_offset=F(1,100))
            c = json.loads((root/'out/certificate.json').read_text())
            r = replay(c)
            self.assertLess(r['width_float'], 1e-7)
            self.assertEqual(r['response_dimension'], 0)
            self.assertEqual(r['reference_action_states'], 3)
            with self.assertRaises(ValueError), contextlib.redirect_stdout(io.StringIO()):
                transfer(source, root/'bad', F(1,5), upper_steps=1, norm_squares=True, temple_offset=F(1,100))
            self.assertFalse((root/'bad/certificate.json').exists())
            self.assertTrue((root/'bad/failure.json').is_file())


if __name__ == '__main__':
    unittest.main()
