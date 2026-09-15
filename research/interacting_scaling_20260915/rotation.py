"""Rational orbital controls and exact two-sided transfer to the original model."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import shutil
import time
from experiments.marginal_symbolic import add, canonical, mono, product, scale
from research.molecular_collective_20260913.core import digest
from research.transfer_followup_20260915.rotated_upper import check_rotation, rational_orthogonal, rotate_and_round
from research.interacting_scaling_20260915.budget import ROOT, OUT, dump


def rotate_polynomial(poly, integer_matrix, denominator):
    check_rotation(integer_matrix, denominator)
    out = {}
    for word, coefficient in canonical(poly).items():
        image = mono((), coefficient)
        for creation, i in word:
            if not 0 <= i < 2*len(integer_matrix):
                raise ValueError('Polynomial does not fit the supplied orbital rotation')
            letter = {((creation, 2*j+i%2),): F(v, denominator)
                for j, v in enumerate(integer_matrix[i//2]) if v}
            image = product(image, letter)
        out = add(out, image)
    return canonical(out)


def construct(source, output):
    import numpy as np
    from pyscf import gto
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    data = json.loads((source/'fixture.json').read_text())
    z = np.load(source/'integrals.npz')
    specification = json.loads((source/'specification.json').read_text()) if (source/'specification.json').exists() else {}
    mol = gto.M(atom=data['geometry'], basis=data['basis'], unit=data['unit'], charge=specification.get('charge', 0), spin=0, verbose=0)
    if mol.nelectron != data['particles']:
        raise ValueError('The orbital model and its declared particle count disagree')
    overlap = mol.intor_symmetric('int1e_ovlp')
    values, vectors = np.linalg.eigh(overlap)
    if min(values) <= 1e-8:
        raise ValueError('Ill-conditioned AO overlap for the declared localization')
    target = z['mo'].T@((vectors*np.sqrt(values))@vectors.T)
    Z, denominator = rational_orthogonal(target)
    rotation = {'integer_matrix': Z, 'denominator': str(denominator),
        'original_fixture_sha256': digest(data), 'design': 'Symmetric AO orthogonalization in atom order'}
    fixture, allowance = rotate_and_round(data, Z, denominator)
    dump(output/'rotation.json', rotation)
    dump(output/'fixture.json', fixture)
    dump(output/'original_fixture.json', data)
    dump(output/'construction.json', {'seconds': time.monotonic()-start,
        'source': str(source), 'source_integrals_sha256': hashlib.sha256((source/'integrals.npz').read_bytes()).hexdigest(),
        'allowance_Ha': str(allowance), 'full_fixed_N_determinants_enumerated': 0})


def initialize_local(name, rotation_directory, widths, complete=False, collective_widths=()):
    from research.interacting_scaling_20260915.dictionary import windows
    case = OUT/'cases'/name
    case.mkdir(parents=True, exist_ok=False)
    data = json.loads((rotation_directory/'fixture.json').read_text())
    for item in ('fixture.json', 'rotation.json'):
        shutil.copyfile(rotation_directory/item, case/item)
    state = rotation_directory/'mps/state.json'
    if state.exists():
        (case/'mps').mkdir()
        shutil.copyfile(state, case/'mps/state.json')
    dump(case/'design.json', {'clusters': [c for width in widths for c in windows(data['modes']//2, width)],
        'collective_clusters': [c for width in collective_widths for c in windows(data['modes']//2, width)],
        'collective_pairs': True, 'complete': complete, 'max_Gram_entries': 2000000,
        'max_coefficient_rows': 180000, 'global_particle_number_only': True,
        'all_declared_Gram_cross_terms_retained': True})
    dump(case/'input_dependencies.json', {'conditional_development_case': True,
        'rotated_input': str(rotation_directory), 'existing_state_reused': state.exists(),
        'source_state_discovery_and_rotation_cost_additional': True,
        'old_Gram_or_coefficient_map_input': False, 'old_winning_factor_input': False})
    return case


def transfer(original, case, accepted, output):
    start = time.monotonic()
    data = json.loads((original/'fixture.json').read_text())
    rotation = json.loads((case/'rotation.json').read_text())
    if rotation['original_fixture_sha256'] != digest(data):
        raise ValueError('Rotation is bound to another original Hamiltonian')
    recomputed, allowance = rotate_and_round(data, rotation['integer_matrix'], int(rotation['denominator']))
    supplied = json.loads((case/'fixture.json').read_text())
    if recomputed != supplied:
        raise ValueError('Transformed Hamiltonian or allowance does not reproduce')
    interval = json.loads((accepted/'interval.json').read_text())
    if json.loads((accepted/'fixture.json').read_text()) != supplied:
        raise ValueError('Accepted lower belongs to another transformed Hamiltonian')
    if not json.loads((accepted/'complete_replay.json').read_text())['all_upper_and_lower_dependencies_rechecked']:
        raise ValueError('Both endpoints must first pass complete exact replay')
    from research.correlated_pair_20260913.mps_exact import check as check_state
    from research.collective_completion_20260914.spin_screen import check as check_lower
    upper_check = check_state(supplied, json.loads((case/'mps/state.json').read_text()))
    lower_check = check_lower(supplied, json.loads((accepted/'certificate.json').read_text()),
        json.loads((accepted/'nonsinglet.json').read_text()))
    if F(upper_check['upper_Ha']) != F(interval['upper_Ha']) or F(lower_check['lower']) != F(interval['lower_Ha']):
        raise ValueError('Transferred endpoints did not reproduce from the actual witnesses')
    L, U = F(interval['lower_Ha'])-allowance, F(interval['upper_Ha'])+allowance
    dump(output, {'lower_Ha': str(L), 'upper_Ha': str(U), 'width_Ha': str(U-L),
        'width_mHa': float(1000*(U-L)), 'target_met': U-L <= F(1, 625),
        'rotation_allowance_each_endpoint_Ha': str(allowance),
        'original_fixture_sha256': digest(data), 'rotated_fixture_sha256': digest(supplied),
        'accepted_interval_sha256': hashlib.sha256((accepted/'interval.json').read_bytes()).hexdigest(),
        'upper_and_lower_witnesses_rechecked': True,
        'argument': 'Exact paired orbital unitary and operator-norm perturbation; each endpoint pays the coefficient allowance.',
        'full_fixed_N_determinants_enumerated': 0, 'seconds': time.monotonic()-start})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='action', required=True)
    c = sub.add_parser('construct'); c.add_argument('source', type=Path); c.add_argument('output', type=Path)
    i = sub.add_parser('initialize'); i.add_argument('name'); i.add_argument('rotation_directory', type=Path)
    i.add_argument('--widths', type=int, nargs='*', default=[2]); i.add_argument('--complete', action='store_true')
    i.add_argument('--collective-widths', type=int, nargs='*', default=[])
    t = sub.add_parser('transfer'); t.add_argument('original', type=Path); t.add_argument('case', type=Path)
    t.add_argument('accepted', type=Path); t.add_argument('output', type=Path)
    a = p.parse_args()
    if a.action == 'construct': construct(a.source.resolve(), a.output.resolve())
    elif a.action == 'initialize': print(initialize_local(a.name, a.rotation_directory.resolve(), a.widths, a.complete, a.collective_widths))
    else: transfer(a.original.resolve(), a.case.resolve(), a.accepted.resolve(), a.output.resolve())
