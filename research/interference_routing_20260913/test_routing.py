from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product as cartesian
from pathlib import Path
import json
import random
import unittest

from experiments.marginal_symbolic import add, canonical, encode, mono, product, scale
from experiments.marginal_hunt_car import adj
from research.all_angles_20260913.frontier.exact_h4_oracle import matrix, ldl_psd, rayleigh
from research.interference_routing_20260913.routing import (
    path_capacity, positive_tree, network, check_network, replay, load_sector,
)
from research.interference_routing_20260913.local_templates import (
    discover, replay as local_replay, RATIOS,
)
from research.interference_routing_20260913.scope_audit import obstruction
from research.interference_routing_20260913.sector_counterexample import EDGES, exchange_matrix, run as sector_counterexample

ROOT = Path(__file__).resolve().parents[2]


def triangle(t):
    return [[1-t, F(-1), t], [F(-1), F(2), F(-1)], [t, F(-1), 1-t]]


def laplacian(n, edges):
    a = [[F(0)]*n for _ in range(n)]
    for i, j, w in edges:
        a[i][i] += w
        a[j][j] += w
        a[i][j] -= w
        a[j][i] -= w
    return a


def car_model(model):
    """Independent CAR expansion; no routing/path logic enters assembly."""
    m = model['modes']
    density = [mono(((1, i), (0, i))) for i in range(m)]
    def swap(i, j):
        parity = mono(())
        for k in range(i+1, j):
            parity = product(parity, add(mono(()), scale(density[k], -2)))
        hop = product(product(mono(((1, i),)), parity), mono(((0, j),)))
        return add(density[i], density[j], scale(product(density[i], density[j]), -2),
                   scale(add(hop, canonical(adj(hop))), -1))
    h = mono((), F(model.get('offset', '0')))
    for i, j in enumerate(model['nearest']):
        h = add(h, scale(swap(i, i+1), F(j)))
    for i, t in enumerate(model['next_nearest']):
        h = add(h, scale(swap(i, i+2), -F(t)))
    for i, field in enumerate(model['fields']):
        h = add(h, scale(density[i], F(field)))
    for i, value in enumerate(model['interaction']):
        h = add(h, scale(product(density[i], density[i+1]), F(value)))
    return {'modes': m, 'hamiltonian': encode(h)}


class RoutingTests(unittest.TestCase):
    def test_frustrated_triangle_and_threshold(self):
        for t in (F(1, 3), F(1, 2)):
            a = triangle(t)
            edges, _ = network(a, [1, 1, 1])
            lower, _ = check_network(a, [1, 1, 1], positive_tree(3, edges))
            self.assertEqual(lower, 0)
            self.assertGreater(a[0][1]*a[1][2]*a[0][2], 0)
        a = triangle(F(3, 5))
        with self.assertRaisesRegex(ValueError, 'not PSD'):
            check_network(a, [1, 1, 1], [[0, 1], [1, 2]])

    def test_shared_routes_are_not_independent(self):
        positive = [(0, 1, F(1)), (1, 2, F(1)), (2, 3, F(1))]
        negative = [(0, 2, F(-1, 2)), (1, 3, F(-1, 2))]
        tree = [[0, 1], [1, 2], [2, 3]]
        for edge in negative:
            k, _ = path_capacity(4, positive+[edge], tree)
            self.assertTrue(ldl_psd(k)[0])
        k, _ = path_capacity(4, positive+negative, tree)
        self.assertFalse(ldl_psd(k)[0])
        a = laplacian(4, positive+negative)
        self.assertEqual(rayleigh(a, [2, 1, -1, -2]), F(-3, 10))

    def test_positive_and_mixed_diagonal_slack(self):
        a = triangle(F(3, 5))
        for i in range(3):
            a[i][i] += F(1, 5)
        with self.assertRaises(ValueError):
            check_network(a, [1, 1, 1], [[0, 1], [1, 2]])
        edges, _ = network(a, [1, 1, 1], F(0))
        self.assertEqual(check_network(a, [1, 1, 1], positive_tree(4, edges), F(0))[0], 0)
        a = [[F(1), F(-2)], [F(-2), F(4)]]
        edges, local = network(a, [1, 1], F(0))
        self.assertEqual(local, [F(-1), F(2)])
        self.assertEqual(check_network(a, [1, 1], positive_tree(3, edges), F(0))[0], 0)
        self.assertEqual(rayleigh(a, [2, 1]), 0)

    def test_congruence_needs_transformed_metric(self):
        self.assertTrue(ldl_psd([[F(2)]])[0])  # THT-bI
        self.assertFalse(ldl_psd([[F(-4)]])[0])  # T(H-bI)T
        self.assertFalse(ldl_psd([[F(0), F(1)], [F(1), F(0)]])[0])

    def test_capacity_reconstructs_every_entry(self):
        rng = random.Random(913)
        for size in range(2, 7):
            tree = [[i, i+1] for i in range(size-1)]
            edges = [(i, j, F(1) if j == i+1 else F(rng.randrange(-3, 4), 5))
                     for i in range(size) for j in range(i+1, size)]
            k, _ = path_capacity(size, edges, tree)
            b = [[int(i == left)-int(i == right) for left, right in tree] for i in range(size)]
            rebuilt = [[sum(b[i][r]*k[r][s]*b[j][s] for r in range(size-1) for s in range(size-1))
                        for j in range(size)] for i in range(size)]
            self.assertEqual(rebuilt, laplacian(size, edges))

    def test_invalid_routes_and_zero_amplitude(self):
        a = triangle(F(1, 3))
        with self.assertRaisesRegex(ValueError, 'nonzero integer'):
            network(a, [1, 0, 1])
        edges, _ = network(a, [1, 1, 1])
        for bad in ([[0, 1]], [[0, 1], [0, 1]], [[0, 1], [0, 2]], [[True, 1], [1, 2]]):
            with self.assertRaises(ValueError):
                path_capacity(3, edges, bad)
        with self.assertRaisesRegex(ValueError, 'Disconnected tree'):
            path_capacity(4, [(0, 1, 1), (1, 2, 1), (0, 2, 1)], [[0, 1], [1, 2], [0, 2]])

    def test_saved_molecular_binding_and_refusals(self):
        fixture = ROOT/'results/certificate_scaling/active_space_ladder/h4/fixture.json'
        proof = json.loads((ROOT/'results/interference_routing_20260913/molecular/h4_line_grounded_proof.json').read_text())
        self.assertLess(F(replay(fixture, proof)['width']), F(1, 625))
        for mutate in (
            lambda p: p['components'].pop(),
            lambda p: p.__setitem__('fixture_sha256', 'bad'),
            lambda p: p.__setitem__('particles', 3),
            lambda p: p['claim'].__setitem__('upper', '-100'),
            lambda p: p['components'][0].__setitem__('lower', '100'),
            lambda p: p.__setitem__('upper_amplitudes', [0]*70),
        ):
            bad = deepcopy(proof)
            mutate(bad)
            with self.assertRaises(ValueError):
                replay(fixture, bad)
        with self.assertRaisesRegex(ValueError, 'Fixture/sector mismatch'):
            replay(ROOT/'results/marginal_molecule/h4_rectangle_sto3g.json', proof)
        with self.assertRaisesRegex(ValueError, 'Enumeration budget'):
            load_sector(ROOT/'results/certificate_scaling/active_space_ladder/h6/fixture.json')


class LocalTemplateTests(unittest.TestCase):
    def model(self, m=6):
        return {'modes': m, 'particles': m//2,
                'nearest': ['3/10' if i%2 == 0 else '3/5' for i in range(m-1)],
                'next_nearest': ['1/10']*(m-2),
                'fields': [str(F(i%3-1, 10)) for i in range(m)],
                'interaction': ['1/7' if i%2 else '-1/11' for i in range(m-1)], 'offset': '2/13'}

    def test_car_all_sectors_and_exact_variational_upper(self):
        model = self.model()
        expanded = car_model(model)
        for particles in range(7):
            model['particles'] = particles
            ratios, _ = discover(model)
            receipt = local_replay({'model': model, 'ratios': ratios})
            a, states = matrix(expanded, particles)
            lower = F(receipt['lower'])
            shifted = [[a[i][j]-(lower if i == j else 0) for j in range(len(a))] for i in range(len(a))]
            self.assertTrue(ldl_psd(shifted)[0], particles)
            self.assertEqual(F(receipt['upper']), rayleigh(a, [1]*len(a)))
            if particles == 3:
                self.assertTrue(any(a[i][j]*a[j][k]*a[i][k] > 0 for i, j, k in combinations(range(len(a)), 3)))

    def test_capacity_ablation_and_tamper(self):
        model = self.model()
        ratios, _ = discover(model)
        result = local_replay({'model': model, 'ratios': ratios})
        with self.assertRaisesRegex(ValueError, 'overdrawn'):
            local_replay({'model': model, 'ratios': ['1']*4})
        for bad_ratios in (['0']*4, ['1']*3, [0.5]*4):
            with self.assertRaises(ValueError):
                local_replay({'model': model, 'ratios': bad_ratios})
        with self.assertRaisesRegex(ValueError, 'False stored'):
            local_replay({'model': model, 'ratios': ratios, 'claim': {'lower': '-100', 'upper': '0', 'width': '100'}})
        self.assertEqual(result['many_body_states_enumerated'], 0)

    def test_support_obstruction_is_not_an_acceptance_rule(self):
        from experiments.marginal_symbolic import decode
        model = self.model()
        h = decode(car_model(model)['hamiltonian'], 6, 4)
        self.assertEqual(obstruction(h), {'status': 'no_support_obstruction', 'is_template_acceptance': False})
        h = add(h, mono(((1, 0), (0, 5)), F(1, 13)))
        result = obstruction(h)
        self.assertEqual(result['status'], 'outside_template_span')
        self.assertEqual(result['witness_coefficient'], '1/13')

    def test_finite_grammar_dp_against_exhaustive_routes(self):
        rng = random.Random(19)
        for _ in range(12):
            model = self.model(5)
            model['nearest'] = [str(F(rng.randrange(1, 9), 10)) for _ in range(4)]
            feasible = False
            for ratios in cartesian(RATIOS, repeat=3):
                used = [F(0)]*4
                for i, rho in enumerate(ratios):
                    used[i] += F(1, 10)*(1+rho)
                    used[i+1] += F(1, 10)*(1+1/rho)
                if all(x <= F(y) for x, y in zip(used, model['nearest'])):
                    feasible = True
                    break
            try:
                ratios, _ = discover(model)
                local_replay({'model': model, 'ratios': ratios})
                found = True
            except ValueError:
                found = False
            self.assertEqual(found, feasible)

    def test_overdrawn_model_has_exact_negative_state(self):
        model = self.model(4)
        model.update(nearest=['1']*3, next_nearest=['1/2']*2,
                     fields=['0']*4, interaction=['0']*3, offset='0', particles=1)
        a, _ = matrix(car_model(model), 1)
        self.assertEqual(rayleigh(a, [2, 1, -1, -2]), F(-3, 10))
        with self.assertRaises(ValueError):
            discover(model)

    def test_extensive_width_is_retained(self):
        model = self.model(64)
        model.update(fields=['0']*64, interaction=['1/10000']*63, offset='0')
        ratios, _ = discover(model)
        result = local_replay({'model': model, 'ratios': ratios})
        self.assertEqual(F(result['width']), F(31, 20000))
        self.assertEqual(result['routing_templates'], 62)

    def test_one_particle_routing_does_not_lift_to_all_particles(self):
        self.assertEqual(sector_counterexample()['two_particle_eigenvalue'], '-2')
        # Independently expand the parity-dressed CAR operator, including the
        # three-link string, and compare every fixed-N matrix entry.
        densities = [mono(((1, i), (0, i))) for i in range(4)]
        h = {}
        for i, j, weight in EDGES:
            parity = mono(())
            for k in range(i+1, j):
                parity = product(parity, add(mono(()), scale(densities[k], -2)))
            hop = product(product(mono(((1, i),)), parity), mono(((0, j),)))
            k_op = add(densities[i], densities[j], scale(product(densities[i], densities[j]), -2),
                       scale(add(hop, canonical(adj(hop))), -1))
            h = add(h, scale(k_op, weight))
        for particles in range(5):
            actual, states = matrix({'modes': 4, 'hamiltonian': encode(h)}, particles)
            expected, expected_states = exchange_matrix(4, particles, EDGES)
            self.assertEqual(states, expected_states)
            self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
