"""Stdlib exact ceiling for two new coefficients with every old field frozen."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import importlib
import json
import sys
from functools import lru_cache

from spin_coherence_obstruction import ROOT, BASE, OLD_MODULES, pure_coherence, read_accepted, hashes
from full_overlap_telescope import moment
from experiments.marginal_local_hubbard_block import _actions, _sector
from experiments.marginal_projector_extendibility import _projector_vector
from experiments.marginal_joint_family_limit import _physical_vector, _weight
from experiments.marginal_charged_projectors import charged_vectors
from experiments.marginal_range_two_density import diagonal_value
from experiments.marginal_quadratic_charge_telescope import local_value as quadratic_value
from experiments.marginal_charge_square_pairs import local_value as square_value
from experiments.marginal_charge_indicator_telescope import local_value as indicator_value
from experiments.marginal_signed_charge_telescope import local_value as signed_value
from experiments.marginal_spin_word_telescope import local_value as word_value

VECTORS = ({346: 1, 409: 1}, {314: 1, 614: 1})


def fixed_expectation(c):
    if c.get('kind') != 'hubbard_projector_extension_v18' or any(F(c['target'][k]) != v for k, v in [('U', 4), ('t', 1), ('V', F(1, 2))]):
        raise ValueError('Matched accepted ENERGYv18 recipe required')
    local = c['local_window']
    u, t, d = [list(map(F, local[key])) for key in ('onsite_profile', 'hopping_profile', 'density_profile')]
    physical = _actions(6, F(10, 3), 1, u, t, F(1, 2), d)
    combined = [dict(image) for image in physical]
    for name in OLD_MODULES:
        module = importlib.import_module('experiments.marginal_' + name)
        field = 'coherent_projector' if name == 'coherent_projector_telescope' else name
        action = module.actions({label: F(value) for label, value in c[field].items()}) if hasattr(module, 'LABELS') else module.actions()
        factor = F(1) if hasattr(module, 'LABELS') else F(c[field])
        for s, image in enumerate(action):
            for target, value in image.items():
                combined[s][target] = combined[s].get(target, F(0)) + factor * value
    sparse = {int(s): F(value) for s, value in c['telescoping_diagonal'].items()}
    range2 = list(map(F, local['range_two_density_profile']))
    fields = [(quadratic_value, c['quadratic_charge_telescope']), (square_value, c['charge_square_pair_telescope']),
              (indicator_value, c['higher_charge_indicator_telescope']), (signed_value, c['signed_charge_telescope']),
              (word_value, c['spin_word_telescope'])]
    fields = [(fn, {key: F(value) for key, value in terms.items()}) for fn, terms in fields]

    @lru_cache(maxsize=4096)
    def diagonal(state):
        return (diagonal_value(state, range2) + sparse.get(state & 1023, F(0)) - sparse.get(state >> 2, F(0))
                + sum((fn(state, terms) for fn, terms in fields), F(0)))

    half, half_norm = _projector_vector(c['vector'], 6)
    charged, charged_norm = charged_vectors(c['joint']['vector'])
    alpha, beta, ratio = F(c['penalty']), F(c['joint']['penalty']), F(c['joint']['ratio'])

    def expectation(vector, norm):
        value = sum((F(a)*b*vector.get(target, 0) for s, a in vector.items()
                     for target, b in combined[s].items()), F(0)) / norm
        value += sum((a*a*diagonal(s) for s, a in vector.items()), F(0))/norm
        ph = F(sum(a*half.get(s, 0) for s, a in vector.items())**2, norm*half_norm)
        q = ph + ratio*sum((F(sum(a*v.get(s, 0) for s, a in vector.items())**2, norm*charged_norm) for v in charged), F(0))
        return value + alpha*ph + beta*q

    return expectation


def replay(c, proposal):
    if proposal.get('kind') != 'spin_coherence_fixed_family_proposal_v1':
        raise ValueError('Explicit frozen-recipe two-coefficient family required')
    mixture = proposal.get('mixture')
    if type(mixture) is not list or not 1 <= len(mixture) <= 3:
        raise ValueError('One through3 physical sources required')
    states = []
    for item in mixture:
        weight = _weight(item['weight'])
        vector, norm = _physical_vector(item['vector'])
        if weight <= 0 or len({_sector(s, 6) for s in vector}) != 1:
            raise ValueError('Positive physical source in one spin sector required')
        states.append((weight, vector, norm))
    if sum(weight for weight, vector, norm in states) != 1:
        raise ValueError('Exact normalized mixture required')
    moments = []
    for vector in VECTORS:
        denominator, _, matrix, _ = pure_coherence(vector)
        value = moment(mixture, matrix, denominator)
        if value:
            raise ValueError('Both new coherence expectations must cancel exactly')
        moments.append(str(value))
    expectation = fixed_expectation(c)
    energies = [expectation(vector, norm) for weight, vector, norm in states]
    average = sum((weight*energy for (weight, vector, norm), energy in zip(states, energies)), F(0))
    offset = F(c['penalty'])*F(c['projector_sum_ceiling'])/c['windows']
    offset += F(c['joint']['penalty'])*F(c['joint']['projector_sum_ceiling'])/c['joint']['windows']
    cap = (average-offset)/5
    lower = (F(c['penalized_lower'])-offset)/5
    if cap < lower:
        raise ValueError('Family cap contradicts accepted seed lower')
    return {'accepted': True, 'mixture_sources': len(states), 'mixture_trace': '1',
            'new_coherence_moments': moments, 'exact_source_local_expectations': list(map(str, energies)),
            'fixed_recipe_family_ceiling': str(cap), 'ceiling_float': float(cap),
            'accepted_seed_lower': str(lower), 'maximum_gain_over_seed': str(cap-lower),
            'maximum_gain_over_seed_float': float(cap-lower),
            'all_old_fields_frozen': True, 'new_coefficients_unrestricted_real': True,
            'proof': 'For every real pair gamma, lambda_min(A+gamma0*T0+gamma1*T1) is at most the expectation in this positive trace-one mixture. Both T expectations vanish exactly, leaving the fixed A expectation. Subtract the fixed projector penalty offset and divide by5. No fidelity ceiling is imposed on this witness because both penalties remain fixed.',
            'scope': 'An exact ceiling only for the two new coefficients with all ENERGYv18 fields frozen. This does not cap reoptimization of older fields, the full enlarged family, or physical ground energy.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    directory = parser.parse_args().directory.resolve()
    proposal_path = directory / 'family_proposal.json'
    proposal = json.loads(proposal_path.read_text())
    cp = ROOT / proposal['seed_certificate']
    files = {proposal_path, cp, Path(__file__).resolve()}
    energy = read_accepted(cp.parent / 'range_two_replay.json', files)
    if str(cp.relative_to(ROOT)) not in energy['source_sha256']:
        raise ValueError('Seed is not bound to accepted energy receipt')
    for case, vector in zip(('W_zero', 'W_plus_1'), VECTORS):
        receipt = read_accepted(BASE / case / 'telescope_replay.json', files)
        if {int(s): a for s, a in receipt['five_site_vector'].items()} != vector:
            raise ValueError('Canonical new coherence source differs')
    result = replay(json.loads(cp.read_text()), proposal)
    if F(result['accepted_seed_lower']) != F(energy['lower_replay']['periodic_lower_density']):
        raise ValueError('Seed lower differs from accepted energy')
    result['source_sha256'] = hashes(files)
    output = directory / 'fixed_family_replay.json'
    if output.exists():
        raise ValueError('Refusing to overwrite exact family receipt')
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('accepted', 'mixture_sources', 'ceiling_float', 'maximum_gain_over_seed_float')}))


if __name__ == '__main__':
    main()
