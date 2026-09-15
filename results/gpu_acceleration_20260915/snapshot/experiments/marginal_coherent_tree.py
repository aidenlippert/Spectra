"""Coherent transition polynomials and exact charge-metric branch coverage.

No determinant action or sector matrix is used. Completed occupation leaves
are counted explicitly: avoiding a matrix is not avoiding enumeration.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product as cartesian
import json
from math import comb
from pathlib import Path

from experiments.marginal_charge_product import feature_orbits
from experiments.marginal_spin_reduction import SpinZeroOracle
from experiments.marginal_symbolic import mono, product


class CoherentCharge:
    def __init__(self, data):
        self.oracle = SpinZeroOracle(data)
        self.modes, self.particles = self.oracle.modes, self.oracle.particles
        if self.modes % 4 or self.particles != self.modes//2:
            raise ValueError('Complete half-filled even-site valence complement required')
        self.sites = self.modes//2; self.all_bits = (1 << self.modes)-1
        self.spin_masks = [sum(1 << i for i in range(spin, self.modes, 2)) for spin in (0, 1)]
        self.target = self.particles//2
        rule = data['metric_rule']
        orbits = feature_orbits(self.sites, rule['family'], rule.get('max_pair_distance'))
        factors = rule['factors']
        if type(factors) is not list or len(factors) != len(orbits) or any(type(x) is not str for x in factors):
            raise ValueError('Exact charge factor for each feature required')
        factors = [F(x) for x in factors]
        if any(x <= 0 for x in factors):
            raise ValueError('Positive charge factors required')
        self.local_factors = [(label, factor) for orbit, factor in zip(orbits, factors) for label in orbit]
        self.groups = {}
        # H monomials with the same changed occupations share a bare CAR
        # transition. Spectator number operators commute with that transition.
        for word, coefficient in self.oracle.h.items():
            creates = {i for c, i in word if c}; annihilates = {i for c, i in word if not c}
            cset, aset = creates-annihilates, annihilates-creates
            if not cset and not aset:
                continue
            if len(cset) != len(aset):
                raise ValueError('Particle-preserving transitions required')
            spectator = creates & annihilates
            if len(spectator) > 1:
                raise ValueError('At most affine spectator amplitudes supported')
            bare = tuple((1, i) for i in sorted(cset))+tuple((0, i) for i in sorted(aset, reverse=True))
            density = {(): F(1)}
            for i in sorted(spectator):
                density = product(density, mono(((1, i), (0, i))))
            expanded = product(mono(bare), density)
            if len(expanded) != 1 or word not in expanded or abs(expanded[word]) != 1:
                raise ValueError('Failed exact CAR spectator factorization')
            c = sum(1 << i for i in cset); a = sum(1 << i for i in aset)
            key = (c, a); mask = sum(1 << i for i in spectator)
            poly = self.groups.setdefault(key, {})
            poly[mask] = poly.get(mask, F(0))+coefficient/expanded[word]
        self.groups = {key: {m: x for m, x in poly.items() if x} for key, poly in self.groups.items()}
        self.groups = {key: poly for key, poly in self.groups.items() if poly}
        self.source_polynomial = data.get('source_polynomial_bound', False)
        if type(self.source_polynomial) is not bool:
            raise ValueError('Explicit Boolean source-polynomial selector required')
        self.joint_line = data.get('joint_line_bound', False)
        if type(self.joint_line) is not bool:
            raise ValueError('Explicit Boolean joint-line selector required')
        self.joint_stats = {'joint_line_branches': 0, 'joint_line_assignment_capacity': 0,
                            'joint_line_ionic_assignments': 0, 'joint_line_metric_endpoint_evaluations': 0,
                            'joint_line_amplitude_endpoint_evaluations': 0,
                            'joint_line_transition_groups_used': 0}
        self.joint_simplex = data.get('joint_simplex_bound', False)
        if type(self.joint_simplex) is not bool:
            raise ValueError('Explicit Boolean joint-simplex selector required')
        self.simplex_stats = {name.replace('joint_line_', 'joint_simplex_'): 0
                              for name in self.joint_stats}
        self.simplex_stats['joint_simplex_maximum_assignments'] = 0
        self.conditioned_sources = set()
        self.max_polynomial_terms = 0
        self.charge_model = None
        if 'diagonal_charge_lambda' in data:
            if type(data['diagonal_charge_lambda']) is not str:
                raise ValueError('Exact charge quadratic threshold required')
            self.charge_model = self.diagonal_charge_model(F(data['diagonal_charge_lambda']))

    def diagonal_charge_model(self, threshold):
        """Reconstruct V=c+ell.q+q^T B q+sum S_ij z_i z_j exactly."""
        from experiments.marginal_charge_polynomial import add, multiply, scale
        from experiments.marginal_schur_transfer import ldl_pivots
        d = self.oracle.diagonal; m = self.sites
        h = [d.get(1 << (2*i), F(0)) for i in range(m)]
        u = [d.get(3 << (2*i), F(0)) for i in range(m)]
        constant = d.get(0, F(0))+sum(h); ell = [h[i]+u[i]/2 for i in range(m)]
        matrix = [[u[i]/2 if i == j else F(0) for j in range(m)] for i in range(m)]
        q = [{0: F(-1), 1 << (2*i): F(1), 1 << (2*i+1): F(1)} for i in range(m)]
        z = [{1 << (2*i): F(1), 1 << (2*i+1): F(-1)} for i in range(m)]
        pair_terms = []; penalty = F(0)
        for i in range(m):
            for j in range(i+1, m):
                same = d.get((1 << (2*i)) | (1 << (2*j)), F(0))
                opposite = d.get((1 << (2*i)) | (1 << (2*j+1)), F(0))
                charge, spin = (same+opposite)/2, (same-opposite)/2
                constant += charge; ell[i] += charge; ell[j] += charge
                matrix[i][j] = matrix[j][i] = charge/2
                matrix[i][i] += abs(spin)/2; matrix[j][j] += abs(spin)/2
                penalty += abs(spin)
                pair_terms += [scale(multiply(q[i], q[j]), charge), scale(multiply(z[i], z[j]), spin)]
        reconstructed = add({0: constant}, *(scale(q[i], ell[i]) for i in range(m)),
                            *(scale(multiply(q[i], q[i]), u[i]/2) for i in range(m)), *pair_terms)
        if reconstructed != {mask: value for mask, value in d.items() if value}:
            raise ValueError('Diagonal charge/spin identity does not match H')
        remainder = [[matrix[i][j]-threshold*(i == j) for j in range(m)] for i in range(m)]
        if ldl_pivots(remainder) is None:
            raise ValueError('Charge quadratic threshold lacks exact positivity')
        return constant-penalty, ell, remainder, threshold

    @lru_cache(maxsize=10000)
    def charge_diagonal_lower(self, mask, bits):
        constant, ell, remainder, threshold = self.charge_model
        choices = [[value for value in range(4) if value & ((mask >> (2*i)) & 3) == ((bits >> (2*i)) & 3)]
                   for i in range(self.sites)]
        anchor = [values[0].bit_count()-1 if len(values) == 1 else 0 for values in choices]
        bounds = []
        for vector in ([0]*self.sites, anchor):
            rv = [sum((row[j]*vector[j] for j in range(self.sites)), F(0)) for row in remainder]
            offset = constant-sum((vector[i]*rv[i] for i in range(self.sites)), F(0))
            linear = [ell[i]+2*rv[i] for i in range(self.sites)]
            dp = {(0, 0, False): F(0)}
            for i, values in enumerate(choices):
                nxt = {}
                for (alpha, beta, ionic), cost in dp.items():
                    for value in values:
                        a, b = alpha+(value & 1), beta+((value >> 1) & 1)
                        if a > self.target or b > self.target:
                            continue
                        charge = value.bit_count()-1
                        key = (a, b, ionic or value == 3)
                        candidate = cost+linear[i]*charge+threshold*charge**2
                        if key not in nxt or candidate < nxt[key]:
                            nxt[key] = candidate
                dp = nxt
            bound = dp.get((self.target, self.target, True))
            if bound is None:
                raise ValueError('Charge lower bound requested on empty ionic branch')
            bounds.append(offset+bound)
        return max(bounds)

    def close(self, mask, bits):
        """Propagate exact fixed spin populations; None means infeasible."""
        for spin_mask in self.spin_masks:
            free = spin_mask & ~mask
            need = self.target-(bits & spin_mask).bit_count()
            if not 0 <= need <= free.bit_count():
                return None
            if need == 0:
                mask |= free
            elif need == free.bit_count():
                mask |= free; bits |= free
        return mask, bits

    def counts(self, mask, bits):
        total = 1
        for spin_mask in self.spin_masks:
            free = (spin_mask & ~mask).bit_count(); need = self.target-(bits & spin_mask).bit_count()
            if not 0 <= need <= free:
                return 0, 0
            total *= comb(free, need)
        # Count valence configurations by site choices and total alpha number.
        dp = {0: 1}
        for site in range(self.sites):
            local_mask = (mask >> (2*site)) & 3; local_bits = (bits >> (2*site)) & 3
            options = [x for x in (1, 2) if x & local_mask == local_bits]
            nxt = {}
            for n, count in dp.items():
                for value in options:
                    k = n+(value == 1); nxt[k] = nxt.get(k, 0)+count
            dp = nxt
        valence = dp.get(self.target, 0)
        return total, total-valence

    def amplitude_range(self, poly, mask, bits):
        constant = poly.get(0, F(0))
        low = high = constant
        for spin_mask in self.spin_masks:
            values = []
            for i in range(self.modes):
                bit = 1 << i
                if not bit & spin_mask:
                    continue
                value = poly.get(bit, F(0))
                if bit & mask:
                    if bit & bits:
                        low += value; high += value
                else:
                    values.append(value)
            need = self.target-(bits & spin_mask).bit_count()
            values.sort()
            low += sum(values[:need], F(0))
            high += sum(values[len(values)-need:], F(0)) if need else F(0)
        return low, high

    @lru_cache(maxsize=100000)
    def ratio_upper(self, delta, mask, bits):
        """Local factor intervals; exact when source site charges are fixed."""
        ratio = F(1)
        for label, factor in self.local_factors:
            if not any(delta[i] for i, _ in label):
                continue
            options = []
            for i, _ in label:
                local_mask = (mask >> (2*i)) & 3; local_bits = (bits >> (2*i)) & 3
                options.append(sorted({value.bit_count()-1 for value in range(4)
                                       if value & local_mask == local_bits}))
            powers = []
            for charges in cartesian(*options):
                before = after = 1
                for charge, (i, power) in zip(charges, label):
                    before *= charge**power; after *= (charge+delta[i])**power
                powers.append(after-before)
            exponent = max(powers) if factor >= 1 else min(powers)
            ratio *= factor**exponent
        return ratio

    def branch_lower(self, mask, bits):
        from experiments.marginal_charge_polynomial import add, indicator, scale
        if self.joint_simplex:
            from experiments.marginal_occupation_simplex import compile_simplex
            compiled = compile_simplex(self, mask, bits)
            if compiled is not None:
                rows, counts = compiled
                self.simplex_stats['joint_simplex_branches'] += 1
                for target, source in [('assignment_capacity','feasible_assignments'), ('ionic_assignments','ionic_assignments'),
                                       ('metric_endpoint_evaluations','metric_scalar_endpoint_evaluations'),
                                       ('amplitude_endpoint_evaluations','amplitude_scalar_endpoint_evaluations'),
                                       ('transition_groups_used','transition_groups_used')]:
                    self.simplex_stats['joint_simplex_'+target] += counts[source]
                self.simplex_stats['joint_simplex_maximum_assignments'] = max(
                    self.simplex_stats['joint_simplex_maximum_assignments'], counts['feasible_assignments'])
                return min(rows)
        if self.joint_line:
            from experiments.marginal_boolean_line import compile_line
            compiled = compile_line(self, mask, bits)
            if compiled is not None:
                row, counts = compiled
                self.joint_stats['joint_line_branches'] += 1
                for target, source in [('assignment_capacity','feasible_assignments'), ('ionic_assignments','ionic_assignments'),
                                       ('metric_endpoint_evaluations','metric_scalar_endpoint_evaluations'),
                                       ('amplitude_endpoint_evaluations','amplitude_scalar_endpoint_evaluations'),
                                       ('transition_groups_used','transition_groups_used')]:
                    self.joint_stats['joint_line_'+target] += counts[source]
                return row[0]+min(F(0), row[1])
        lower = self.oracle.diagonal_lower(mask, bits)
        if self.charge_model is not None:
            lower = max(lower, self.charge_diagonal_lower(mask, bits))
        polynomial = {}
        if self.source_polynomial:
            for support, value in self.oracle.diagonal.items():
                if not support & mask & ~bits:
                    polynomial = add(polynomial, {support & ~mask: value})
        for (c, a), poly in self.groups.items():
            required = c | a
            if (bits ^ a) & mask & required:
                continue
            closed = self.close(mask | required, bits | a)
            if closed is None:
                continue
            fixed, source = closed
            _, sources = self.counts(fixed, source)
            if not sources:
                continue
            if fixed == self.all_bits:
                self.conditioned_sources.add(source)
            # The transition is a bijection on this source branch. If its
            # entire image is valence, it contributes nothing to QHQ.
            _, targets = self.counts(fixed, source ^ required)
            if not targets:
                continue
            lo, hi = self.amplitude_range(poly, fixed, source)
            magnitude = max(abs(lo), abs(hi))
            if not magnitude:
                continue
            delta = tuple(((c >> (2*i)) & 3).bit_count()-((a >> (2*i)) & 3).bit_count()
                          for i in range(self.sites))
            penalty = magnitude*self.ratio_upper(delta, fixed, source)
            lower -= penalty
            if self.source_polynomial:
                polynomial = add(polynomial, scale(indicator(required & ~mask, a & ~mask), -penalty))
        if self.source_polynomial:
            self.max_polynomial_terms = max(self.max_polynomial_terms, len(polynomial))
            lower = max(lower, self.cardinality_polynomial_lower(polynomial, mask, bits))
        return lower

    def cardinality_polynomial_lower(self, polynomial, mask, bits):
        """For each spin bidegree, sum its smallest required coefficients.

Exactly C(k_alpha,a)C(k_beta,b) monomials of bidegree(a,b) equal one.
Missing coefficients are zeros. No occupation assignments are enumerated.
"""
        free = [(s & ~mask).bit_count() for s in self.spin_masks]
        need = [self.target-(s & bits).bit_count() for s in self.spin_masks]
        groups = {}
        for support, coefficient in polynomial.items():
            if support & mask:
                raise ValueError('Polynomial must use only free occupation bits')
            degree = tuple((support & s).bit_count() for s in self.spin_masks)
            if coefficient:
                groups.setdefault(degree, []).append(coefficient)
        lower = F(0)
        for (a, b), values in groups.items():
            count = (comb(need[0], a) if a <= need[0] else 0)*(comb(need[1], b) if b <= need[1] else 0)
            if not count:
                continue
            total = comb(free[0], a)*comb(free[1], b)
            negatives = sorted(x for x in values if x < 0)
            positives = sorted(x for x in values if x > 0)
            used = min(count, len(negatives)); lower += sum(negatives[:used], F(0))
            remaining = max(0, count-used-(total-len(values)))
            lower += sum(positives[:remaining], F(0))
        return lower

    def cover(self, gamma, tree=None, max_nodes=10000, *, chart_construction=False):
        gamma = F(gamma); building = tree is None
        if type(chart_construction) is not bool or (chart_construction and
                (not building or not (self.joint_line or self.joint_simplex))):
            raise ValueError('Chart construction requires building with a joint chart bound')
        if type(max_nodes) is not int or not 1 <= max_nodes <= 100000:
            raise ValueError('Bounded positive tree node budget required')
        stats = {'nodes': 0, 'bound_leaves': 0, 'exact_occupation_leaves': 0,
                 'pruned_Q_configurations': 0, 'covered_Q_configurations': 0, 'valence_only_leaves': 0}
        def visit(mask, bits, node):
            stats['nodes'] += 1
            if stats['nodes'] > max_nodes:
                raise ValueError('Coherent branch coverage exhausted node budget')
            closed = self.close(mask, bits)
            if closed is None:
                raise ValueError('Tree contains infeasible spin branch')
            mask, bits = closed
            total, qcount = self.counts(mask, bits)
            if not qcount:
                if not building and node != ['valence']:
                    raise ValueError('Incorrect valence-only leaf')
                stats['valence_only_leaves'] += 1
                return ['valence']
            probe = True
            if building and chart_construction and mask != self.all_bits:
                free = [i for i in range(self.modes) if not mask & (1 << i)]
                probe = False
                if self.joint_simplex:
                    from experiments.marginal_occupation_simplex import chart_modes
                    probe = chart_modes(self, mask, bits) is not None
                if self.joint_line and len(free) == 2 and free[0] % 2 == free[1] % 2:
                    probe = self.target-(bits & self.spin_masks[free[0] % 2]).bit_count() == 1
            if (building and probe) or node == ['bound']:
                lower = self.branch_lower(mask, bits)
                if lower >= gamma:
                    stats['bound_leaves'] += 1; stats['covered_Q_configurations'] += qcount
                    if mask == self.all_bits:
                        stats['exact_occupation_leaves'] += 1
                    else:
                        stats['pruned_Q_configurations'] += qcount
                    return ['bound']
                if not building or mask == self.all_bits:
                    raise ValueError('Coherent occupation bound misses requested threshold')
            if building:
                bit = (self.all_bits ^ mask) & -(self.all_bits ^ mask)
                index = bit.bit_length()-1
                node = ['split', index, None, None]
            if (type(node) is not list or len(node) != 4 or node[0] != 'split'
                    or type(node[1]) is not int or not 0 <= node[1] < self.modes or mask & (1 << node[1])):
                raise ValueError('Invalid complete binary occupation split')
            bit = 1 << node[1]
            return ['split', node[1], visit(mask | bit, bits, node[2]), visit(mask | bit, bits | bit, node[3])]
        actual = visit(0, 0, tree)
        expected = self.oracle.sector_dimension-comb(self.sites, self.target)
        if stats['covered_Q_configurations'] != expected:
            raise ValueError('Incomplete ionic-sector coverage')
        result = dict(stats, complement_lower=str(gamma), determinant_actions=len(self.oracle.cache),
                            coherent_transition_groups=len(self.groups), source_hamiltonian_terms=len(self.oracle.h),
                            diagonal_charge_bound=self.charge_model is not None,
                            source_polynomial_bound=self.source_polynomial,
                            maximum_source_polynomial_terms=self.max_polynomial_terms,
                            fully_conditioned_transition_sources=len(self.conditioned_sources))
        if self.joint_line:
            result.update(joint_line_bound=True, **self.joint_stats)
        if self.joint_simplex:
            result.update(joint_simplex_bound=True, **self.simplex_stats)
        if chart_construction:
            result['construction_chart_only'] = True
        return actual, result


def replay(certificate):
    if certificate.get('kind') != 'valence_coherent_tree_v1':
        raise ValueError('Unsupported coherent charge tree')
    if type(certificate.get('tree')) is not list:
        raise ValueError('Explicit occupation coverage tree required')
    oracle = CoherentCharge(certificate)
    _, result = oracle.cover(F(certificate['target_lower']), certificate['tree'], certificate['max_nodes'])
    result['scope'] = ('Exact grouped-CAR amplitudes, unrounded positive charge-product ratios and complete fixed-spin occupation branch coverage of the ionic complement. No sector matrix or determinant action is used. Exact occupation leaves remain explicit enumeration; no efficient scaling claim.')
    return result


def build(source, output, target, max_nodes=10000, *, charge_diagonal=False, source_polynomial=False,
          joint_line=False, joint_simplex=False, chart_construction=False):
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError('Preserve previous coherent tree export')
    data = json.loads(source.read_text())
    cert = {k: data[k] for k in ('modes', 'particles', 'hamiltonian', 'metric_rule')}
    if source_polynomial:
        cert['source_polynomial_bound'] = True
    if joint_line:
        cert['joint_line_bound'] = True
    if joint_simplex:
        cert['joint_simplex_bound'] = True
    oracle = CoherentCharge(cert)
    if charge_diagonal:
        import numpy as np
        # A sufficiently negative trial threshold obtains the exact matrix;
        # the final short rational threshold is checked independently again.
        trial = -1-sum(abs(value) for value in oracle.oracle.diagonal.values())
        _, _, remainder, _ = oracle.diagonal_charge_model(trial)
        matrix = [[float(remainder[i][j]+trial*(i == j)) for j in range(oracle.sites)] for i in range(oracle.sites)]
        threshold = F(int(np.floor(np.linalg.eigvalsh(matrix)[0]*10**6))-1, 10**6)
        cert['diagonal_charge_lambda'] = str(threshold)
        oracle = CoherentCharge(cert)
    tree, construction = oracle.cover(F(target), max_nodes=max_nodes, chart_construction=chart_construction)
    cert.update(kind='valence_coherent_tree_v1', target_lower=str(F(target)), tree=tree, max_nodes=max_nodes)
    receipt = replay(cert)
    output.mkdir(parents=True)
    for name, value in [('certificate', cert), ('receipt', receipt), ('construction_receipt', construction)]:
        (output/(name+'.json')).write_text(json.dumps(value, indent=2)+'\n')
    return receipt


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('--verify'); p.add_argument('--source'); p.add_argument('--output')
    p.add_argument('--target', default='-6.24'); p.add_argument('--max-nodes', type=int, default=10000)
    p.add_argument('--charge-diagonal', action='store_true')
    p.add_argument('--source-polynomial', action='store_true')
    p.add_argument('--joint-line', action='store_true')
    p.add_argument('--joint-simplex', action='store_true')
    p.add_argument('--chart-construction', action='store_true')
    args = p.parse_args()
    print(json.dumps(replay(json.loads(Path(args.verify).read_text())) if args.verify
                     else build(args.source, args.output, F(args.target), args.max_nodes, charge_diagonal=args.charge_diagonal,
                                source_polynomial=args.source_polynomial, joint_line=args.joint_line,
                                joint_simplex=args.joint_simplex, chart_construction=args.chart_construction), indent=2))
