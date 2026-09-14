"""Exact pure-coherence constraints after complete diagonal spin-word closure."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
import hashlib
import importlib
import json
import sys

from full_overlap_telescope import ROOT, build, embed, conjugate, shift6, moment, digest
from experiments.marginal_local_hubbard_block import _sector, _reflection
from experiments.marginal_charged_projectors import charged_vectors, _particlehole, _spinflip

BASE = ROOT / 'results/marginal_graded_hubbard8/residual_coherence'
SELECTION = {'W_zero': 'refined', 'W_plus_1': 'scaled'}
OLD_MODULES = ('hopping_telescope', 'spin_telescope', 'spectator_hopping',
               'pair_transfer', 'two_spectator_hopping', 'three_spectator_hopping',
               'coherent_projector_telescope', 'pure_coherence')


def pure_coherence(vector):
    denominator, full_y, full_t = build(vector)
    y = {key: value for key, value in full_y.items() if key[0] != key[1]}
    t = {key: value for key, value in full_t.items() if key[0] != key[1]}
    if not y or not t:
        raise ValueError('Nonzero offdiagonal coherence required')
    left, right = embed(y, True), embed(y, False)
    expected = {key: left.get(key, 0) - right.get(key, 0) for key in left.keys() | right.keys()}
    if t != {key: value for key, value in expected.items() if value}:
        raise ValueError('Removing diagonal failed to commute with embedding')
    if conjugate(left, shift6) != right:
        raise ValueError('Fermionic translation failed')
    for (r, c), value in t.items():
        if r == c or t.get((c, r)) != value or _sector(r, 6) != _sector(c, 6):
            raise ValueError('Pure coherence is not Hermitian and spin preserving')
    for action in (_particlehole, _spinflip, lambda s: _reflection(s, 6)):
        if conjugate(t, action) != t:
            raise ValueError('Pure coherence breaks a required symmetry')
    total = defaultdict(int)
    translated = t
    for _ in range(6):
        for key, value in translated.items():
            total[key] += value
        translated = conjugate(translated, shift6)
    if translated != t or any(total.values()):
        raise ValueError('Full periodic telescope cancellation failed')
    rows = defaultdict(int)
    for (r, c), value in y.items():
        rows[r] += abs(value)
    y_bound = F(max(rows.values()), denominator)
    return denominator, y, t, y_bound


def functional(matrix, terms):
    return sum(F(sign) * matrix.get((r, c), 0) for r, c, sign in terms)


def eligible_entry(r, c):
    if not r < c or _sector(r, 6) != (4, 2) or _sector(c, 6) != (4, 2):
        return False
    bits = [i for i in range(12) if (r ^ c) >> i & 1]
    return len(bits) != 2 or bits[1] - bits[0] != 2


def discover_functionals(actions, matrices):
    """Exact row elimination; retain coordinate combinations annihilating old columns."""
    count = len(actions)
    if count != 75 or len(matrices) != 2:
        raise ValueError('Bounded75-old plus2-new system required')
    coordinates = {key for matrix in matrices for key in matrix if eligible_entry(*key)}
    for label, action in actions:
        coordinates.update((r, c) for c, image in enumerate(action) for r in image if eligible_entry(r, c))
    if len(coordinates) > 10000:
        raise ValueError('Coordinate elimination budget exceeded')
    basis, found = {}, []
    for coordinate in sorted(coordinates):
        r, c = coordinate
        row = [F(action[c].get(r, 0)) for _, action in actions]
        row += [F(matrix.get(coordinate, 0)) for matrix in matrices]
        functional = {coordinate: F(1)}
        for pivot, (old, terms) in sorted(basis.items()):
            factor = row[pivot]
            if not factor:
                continue
            row = [a-factor*b for a,b in zip(row, old)]
            for key, value in terms.items():
                functional[key] = functional.get(key, F(0))-factor*value
                if not functional[key]:
                    del functional[key]
        pivot = next((i for i,v in enumerate(row) if v), None)
        if pivot is None:
            continue
        factor = row[pivot]
        row = [v/factor for v in row]
        functional = {k:v/factor for k,v in functional.items()}
        basis[pivot] = row, functional
        if pivot >= count:
            found.append(tuple((r,c,a) for (r,c),a in sorted(functional.items())))
        if len(found) == 2:
            return tuple(found)
    raise ValueError('Two independent new directions were not established')


def verify_annihilators(actions, matrices, functionals):
    for label, action in actions:
        if len(action) != 4096:
            raise ValueError('Complete4096-state old action required')
        for terms in functionals:
            if sum(sign * action[c].get(r, 0) for r, c, sign in terms):
                raise ValueError('Functional does not annihilate old term: ' + label)
    minor = [[functional(matrix, terms) for matrix in matrices] for terms in functionals]
    determinant = minor[0][0] * minor[1][1] - minor[0][1] * minor[1][0]
    if not determinant:
        raise ValueError('New directions are not independent modulo the old family')
    return minor, determinant


def read_accepted(path, files):
    data = json.loads(path.read_text())
    if data.get('accepted') is not True:
        raise ValueError('Accepted source receipt required')
    files.add(path)
    for name, expected in data['source_sha256'].items():
        source = ROOT / name
        if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise ValueError('Stale accepted source: ' + name)
        files.add(source)
    return data


def hashes(files):
    files.add(Path(__file__).resolve())
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.is_relative_to(ROOT / 'experiments') or path.name in ('full_overlap_telescope.py', 'full_overlap_density.py'):
                files.add(path)
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(files)}


def main():
    BASE.mkdir(parents=True, exist_ok=False)
    files, cases, matrices, families = set(), {}, [], {}
    for case, selection in SELECTION.items():
        source = BASE.parent / 'pure_coherence' / case / selection
        cp = source / 'range_two_family_limit_certificate.json'
        family = json.loads(cp.read_text())
        family_receipt = read_accepted(source / 'range_two_family_limit_replay.json', files)
        closure = read_accepted(source / 'pure_coherence_closure.json', files)
        overlap = read_accepted(source / 'overlap/full_overlap_replay.json', files)
        if not family_receipt.get('pure_coherence') or not family_receipt.get('full_spin_word') or not closure.get('full_spin_word_overlap_exactly_zero') or overlap.get('difference_diagonal_nonzeros') != 0:
            raise ValueError('Exact complete diagonal closure required before removing diagonal')
        for receipt in (family_receipt, closure, overlap):
            if str(cp.relative_to(ROOT)) not in receipt['source_sha256']:
                raise ValueError('Unmatched family source')
        vector = {int(s): a for s, a in overlap['witness']['five_site_vector'].items()}
        denominator, y, t, y_bound = pure_coherence(vector)
        original_denominator, _, original = build(vector)
        full_moment = moment(family['mixture'], original, original_denominator)
        pure_moment = moment(family['mixture'], t, denominator)
        if not pure_moment or pure_moment != full_moment or pure_moment != F(overlap['witness']['exact_difference']):
            raise ValueError('Pure coherence does not match independently reconstructed overlap')
        if {_sector(int(s), 6) for s, a in family['half_vector'].items() if a} != {(3, 3)}:
            raise ValueError('Fixed half projector has unexpected sector')
        charged, _ = charged_vectors(family['charged_vector'])
        if any(s.bit_count() not in (5, 7) for v in charged for s in v):
            raise ValueError('Fixed charged projectors have unexpected sectors')
        matrix = {key: F(value, denominator) for key, value in t.items()}
        matrices.append(matrix)
        families[case] = family['mixture']
        cases[case] = {
            'accepted': True, 'energy_certificate': False, 'family_directory': str(source.relative_to(ROOT)),
            'five_site_vector': {str(s): a for s, a in vector.items()}, 'matrix_denominator': denominator,
            'five_site_numerator': [[r, c, a] for (r, c), a in sorted(y.items())],
            'six_site_numerator': [[r, c, a] for (r, c), a in sorted(t.items())],
            'five_site_nonzeros': len(y), 'six_site_nonzeros': len(t),
            'six_site_numerator_sha256': digest(t), 'five_site_norm_bound': str(y_bound),
            'six_site_and_open_boundary_norm_bound': str(2 * y_bound),
            'changed_fermion_bit_counts': sorted({(r ^ c).bit_count() for r, c in t}),
            'exact_moment': str(pure_moment), 'moment_float': float(pure_moment),
            'removed_diagonal_moment': str(full_moment - pure_moment),
            'zero_diagonal': True, 'hermitian_spin_preserving': True,
            'all_required_symmetries_and_fermionic_translation_checked': True,
            'six_site_periodic_sum_exactly_zero': True,
            'norm_proof': 'For Hermitian Y, the maximum absolute row sum bounds its spectral norm. Both embedded copies have that norm. Their difference and any open translated telescoping sum have norm at most twice the Y bound.',
            'scope': 'Exact violated stationary coherence constraint after complete diagonal closure. The diagonal part can be removed for these accepted mixtures. No new energy bound, family ceiling or general representability conclusion.',
        }
    old, counts = [], {}
    for name in OLD_MODULES:
        module = importlib.import_module('experiments.marginal_' + name)
        labels = getattr(module, 'LABELS', ('single',))
        counts[name] = len(labels)
        for label in labels:
            action = module.actions({label: 1}) if hasattr(module, 'LABELS') else module.actions()
            old.append((name + ':' + label, action))
    functionals = discover_functionals(old, matrices)
    for terms in functionals:
        for r, c, sign in terms:
            if not eligible_entry(r, c):
                raise ValueError('Functional cannot exclude diagonal, one-body and fixed projector terms')
    minor, determinant = verify_annihilators(old, matrices, functionals)
    cross = {case: [str(moment(mixture, matrix, 1)) for matrix in matrices] for case, mixture in families.items()}
    common = hashes(files)
    for case, data in cases.items():
        data['source_sha256'] = common
        directory = BASE / case
        directory.mkdir()
        (directory / 'telescope_replay.json').write_text(json.dumps(data, indent=2) + '\n')
    support = {
        'accepted': True, 'prior_telescope_counts': counts, 'prior_telescope_count': len(old),
        'annihilator_functionals': [[[r,c,str(a)] for r,c,a in terms] for terms in functionals], 'new_direction_minor': [[str(v) for v in row] for row in minor],
        'new_direction_minor_determinant': str(determinant), 'independent_new_directions': 2,
        'cross_moments_on_selected_mixtures': cross,
        'support_proof': 'Every functional cancels each of the75 prior nondiagonal actions exactly. Its entries are offdiagonal and exclude every nearest-neighbor one-body transition, so every diagonal and physical one-body term vanishes. Sector(4,2) excludes both actual fixed-projector families. A nonsingular2x2 minor proves two independent new directions modulo the entire ENERGYv19/FAMILYv15 affine family.',
        'source_sha256': common,
        'scope': 'Independence of new pure-coherence directions; no claim of improved energy or an attained numerical optimum.',
    }
    (BASE / 'support_replay.json').write_text(json.dumps(support, indent=2) + '\n')
    print(json.dumps({'accepted': True, 'old_directions': len(old), 'new_directions': 2,
                      'minor_determinant': str(determinant),
                      'cases': {k: {n: v[n] for n in ('five_site_nonzeros', 'six_site_nonzeros', 'moment_float', 'six_site_and_open_boundary_norm_bound')} for k, v in cases.items()}}))


if __name__ == '__main__':
    main()
