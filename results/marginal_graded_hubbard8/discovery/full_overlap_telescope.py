"""Exact symmetry-compatible telescope from a sparse positive five-site projector.

Discovery receipt only: this does not change the energy certificate schema.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

from full_overlap_density import ROOT, canonical, transformed
from experiments.marginal_local_hubbard_block import _reflection, _sector
from experiments.marginal_charged_projectors import _particlehole, _spinflip


def permutation(state, modes):
    occupied = [modes[m] for m in range(len(modes)) if state >> m & 1]
    sign = (-1)**sum(a > b for i, a in enumerate(occupied) for b in occupied[i+1:])
    return sum(1 << m for m in occupied), sign


def spin5(state):
    return permutation(state, [m ^ 1 for m in range(10)])


def hole5(state):
    return 1023 ^ state, (-1)**sum(m+m//2 for m in range(10) if state >> m & 1)


def shift6(state):
    return permutation(state, [(m+2) % 12 for m in range(12)])


def add_outer(matrix, vector, sign):
    for column, b in vector.items():
        for row, a in vector.items():
            matrix[row, column] += sign*a*b


def conjugate(matrix, action):
    result = {}
    for (row, column), value in matrix.items():
        r, a = action(row); c, b = action(column)
        result[r, c] = a*b*value
    return result


def embed(matrix, left):
    result = {}
    for (row, column), value in matrix.items():
        for exterior in range(4):
            r, c = ((row | exterior << 10, column | exterior << 10) if left
                    else (row << 2 | exterior, column << 2 | exterior))
            result[r, c] = value
    return result


def build(vector):
    if (type(vector) is not dict or not 1 <= len(vector) <= 2
            or any(type(s) is not int or not 0 <= s < 1024
                   or type(a) is not int or not 0 < abs(a) <= 10**9 for s, a in vector.items())
            or len({_sector(s, 5) for s in vector}) != 1):
        raise ValueError('One or two bounded integer amplitudes in one five-site spin sector required')
    norm = sum(a*a for a in vector.values())
    images = [vector]
    for action in (hole5, spin5):
        images += [transformed(v, action) for v in list(images)]
    orbit = Counter(canonical(v) for v in images)
    for action in (hole5, spin5):
        actual = Counter()
        for v, multiplicity in orbit.items():
            actual[canonical(transformed(dict(v), action))] += multiplicity
        if actual != orbit:
            raise ValueError('Five-site density orbit is not closed')
    odd = defaultdict(int)
    for v in images:
        if sum(a*a for a in v.values()) != norm:
            raise ValueError('Projector image norm changed')
        add_outer(odd, v, 1)
        add_outer(odd, transformed(v, lambda s: _reflection(s, 5)), -1)
    odd = {key: value for key, value in odd.items() if value}
    denominator = 8*norm
    if conjugate(odd, lambda s: _reflection(s, 5)) != {key: -v for key, v in odd.items()}:
        raise ValueError('Five-site source is not reflection odd')
    left, right = embed(odd, True), embed(odd, False)
    if conjugate(left, shift6) != right:
        raise ValueError('Fermionic translation does not give the right embedding')
    telescope = {key: left.get(key, 0)-right.get(key, 0) for key in left.keys() | right.keys()}
    telescope = {key: value for key, value in telescope.items() if value}
    for (row, column), value in telescope.items():
        if telescope.get((column, row)) != value or _sector(row, 6) != _sector(column, 6):
            raise ValueError('Telescope is not Hermitian or changes spin numbers')
    for action in (_particlehole, _spinflip, lambda s: _reflection(s, 6)):
        if conjugate(telescope, action) != telescope:
            raise ValueError('Six-site telescope breaks a required symmetry')
    # Sum all six fermionic cyclic translates; verify cancellation on the full Fock space.
    summed = defaultdict(int); current = telescope
    for _ in range(6):
        for key, value in current.items():
            summed[key] += value
        current = conjugate(current, shift6)
    if current != telescope or any(summed.values()):
        raise ValueError('Periodic telescope identity failed')
    return denominator, odd, telescope


def moment(mixture, matrix, denominator):
    answer = F(0)
    for item in mixture:
        vector = {int(s): a for s, a in item['vector'].items() if a}
        value = sum(a*vector.get(row, 0)*vector.get(column, 0)
                    for (row, column), a in matrix.items())
        answer += F(item['weight'])*F(value, denominator*sum(a*a for a in vector.values()))
    return answer


def digest(matrix):
    return hashlib.sha256(''.join(f'{r},{c}:{v}\n' for (r, c), v in sorted(matrix.items())).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('family_directory', type=Path)
    parser.add_argument('overlap_directory', type=Path)
    args = parser.parse_args(); started = time.monotonic()
    cp = args.family_directory/'range_two_family_limit_certificate.json'
    rp = args.overlap_directory/'full_overlap_replay.json'
    family = json.loads(cp.read_text()); overlap = json.loads(rp.read_text())
    if not overlap.get('accepted') or not overlap.get('stationary_extension_refuted'):
        raise ValueError('Accepted nonstationary full-overlap receipt required')
    for name, expected in overlap['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Stale full-overlap receipt source')
    if str(cp.resolve().relative_to(ROOT)) not in overlap['source_sha256']:
        raise ValueError('The supplied family is not the overlap source')
    vector = {int(s): a for s, a in overlap['witness']['five_site_vector'].items()}
    denominator, odd, telescope = build(vector)
    expectation = moment(family['mixture'], telescope, denominator)
    if not expectation or expectation != F(overlap['witness']['exact_difference']):
        raise ValueError('Direct symmetry-compatible moment differs from full-overlap witness')
    independent = [(r, c, a) for (r, c), a in sorted(telescope.items())
                   if r < c and _sector(r, 6) == (4, 2) and (r ^ c).bit_count() == 8]
    if not independent:
        raise ValueError('No eight-bit matrix entry in the (4,2) sector establishes new support')
    row, column, value = independent[0]
    files = {Path(__file__).resolve(), cp.resolve(), rp.resolve(),
             Path(sys.modules['full_overlap_density'].__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name and str(Path(name).resolve()).startswith(str(ROOT/'experiments')+'/'):
            files.add(Path(name).resolve())
    result = {
        'accepted': True, 'energy_certificate': False,
        'five_site_vector': {str(s): a for s, a in vector.items()},
        'definition': 'A=(P+HPH*+SPS*+SH P H*S*)/4; Y=(A-RAR*)/2; T=Y_left-Y_right; P=vv*/||v||^2.',
        'matrix_denominator': denominator,
        'five_site_numerator': [[r, c, a] for (r, c), a in sorted(odd.items())],
        'five_site_nonzeros': len(odd), 'six_site_nonzeros': len(telescope),
        'six_site_numerator_sha256': digest(telescope),
        'exact_direct_mixture_expectation': str(expectation), 'expectation_float': float(expectation),
        'matches_full_overlap_expectation': True,
        'all_4096_fock_states_covered': True,
        'hermitian': True, 'spin_number_preserving': True,
        'particle_hole_spin_flip_reflection_invariant': True,
        'fermionic_left_to_right_translation_checked': True,
        'six_site_periodic_sum_exactly_zero': True,
        'operator_norm_bound': 1,
        'norm_proof': 'A and RAR* are averages of normalized positive projectors, hence each is between0 andI. Each difference of its two embeddings has norm at most1. T is half the difference of those two differences, so ||T||<=1. Open translated sums leave two Y boundary embeddings and obey the same bound.',
        'new_support_entry': {'row': row, 'column': column, 'sector': [4, 2],
                              'changed_fermion_bits': 8, 'matrix_element': str(F(value, denominator))},
        'support_comparison': 'Diagonal terms change0 bits, one-hop terms (including every charge-spectator multiplier) change2, spin exchange and pair transfer change4. The fixed half-filled projector is supported in (3,3); charged projectors have total particle number5 or7. Thus every prior affine-family term vanishes at this (4,2), eight-bit entry.',
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)},
        'seconds': time.monotonic()-started,
        'scope': 'Exact compact stationary constraint violated by this particular positive local mixture. Independently evaluated on original mixture sources and matched to full symmetry-averaged RDM witness. No new lower energy bound, family ceiling, convergence or general representability claim.'}
    output = args.overlap_directory/'telescope_replay.json'
    if output.exists():
        raise ValueError('Refusing to overwrite a telescope receipt')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('accepted', 'five_site_nonzeros', 'six_site_nonzeros',
                                                 'expectation_float', 'new_support_entry', 'seconds')}), flush=True)


if __name__ == '__main__':
    main()
