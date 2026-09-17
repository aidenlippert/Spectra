"""Exact constructive redundancy check against the current frozen dictionaries."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
from experiments.marginal_symbolic import add, adj, canonical, encode, mono, product
from research.acceptance_channels_20260915.channel import supplied_channel
from research.acceptance_channels_20260915.campaign import OUT, ROOT, dump
from research.interacting_scaling_20260915.spin_patterns import twirl


def coordinates(polynomial, group):
    dictionary = {}
    for i, raw in enumerate(group['words']):
        p = canonical(mono(tuple(tuple(letter) for letter in raw)))
        if len(p) != 1:
            raise ValueError('Expected a monomial dictionary')
        word, sign = next(iter(p.items()))
        if word in dictionary:
            raise ValueError('Repeated canonical word')
        dictionary[word] = (i, sign)
    if not set(polynomial) <= dictionary.keys():
        return None
    vector = [Fraction(0)] * len(group['words'])
    for word, coefficient in polynomial.items():
        i, sign = dictionary[word]
        vector[i] = coefficient / sign
    # This vector defines an exact rank-one PSD addition v v^T in that block.
    return vector


def run():
    B, _, P = supplied_channel()
    adjoint = adj(B)
    rotated_adjoint = canonical({tuple((c, i ^ 1) for c, i in w): a for w, a in adjoint.items()})
    replacement = add(product(adj(B), B), product(adj(rotated_adjoint), rotated_adjoint))
    if twirl(P) != twirl(replacement):
        raise ValueError('Spin-rotated highest-weight replacement does not reproduce the twirled separator')
    results = []
    for name in ('h8_matched_direct_level0', 'h12_main_refined'):
        prepared = ROOT/'results/interacting_scaling_20260915/cases'/name/'prepared'
        frame = json.loads((prepared/'frame.json').read_text())
        bases = np.load(prepared/'bases.npz')
        witnesses = []
        for label, p in (('B', B), ('spin_rotated_B_adjoint', rotated_adjoint)):
            found = None
            for k, block in enumerate(frame['blocks']):
                group = frame['groups'][block['physical_group']]
                if block.get('kind') != 'spin_three_half_highest':
                    continue
                vector = coordinates(p, group)
                if vector is None:
                    continue
                V = bases[f'V_{k}']
                if not np.array_equal(V, np.eye(len(group['words']))):
                    raise ValueError('The retained highest-weight basis is not exactly the identity')
                recovered = add(*(mono(tuple(tuple(a) for a in raw), v)
                    for raw, v in zip(group['words'], vector) if v))
                if canonical(recovered) != canonical(p):
                    raise ValueError('Explicit dictionary coordinates failed')
                found = {'operator': label, 'block': k, 'group': group['name'],
                    'dimension': block['dimension'],
                    'nonzero_coordinates': {str(i): str(v) for i, v in enumerate(vector) if v},
                    'rank_one_PSD_Gram': 'v v^T with these exact rational coordinates'}
                break
            if found is None:
                raise ValueError((name, label, 'operator is not in a retained identity basis'))
            witnesses.append(found)
        results.append({'case': name, 'frame_sha256': hashlib.sha256((prepared/'frame.json').read_bytes()).hexdigest(),
            'basis_sha256': hashlib.sha256((prepared/'bases.npz').read_bytes()).hexdigest(),
            'witnesses': witnesses, 'twirled_separator_already_in_current_cone': True})
    dump(OUT/'channel_membership.json', {
        'kind': 'constructive_exact_cone_membership',
        'spin_rotated_B_adjoint': encode(rotated_adjoint),
        'twirled_polynomial_identity_checked_exactly': True,
        'old_external_paired_family_is_a_different_family': True,
        'external_dual_not_replayed': True,
        'operator_coordinate_convention': 'The supplied polynomial is interpreted verbatim in each current dictionary orbital labelling.',
        'external_orbital_basis_match_established': False,
        'transport_of_external_physical_operator_checked': False,
        'conclusion': 'Adding this isolated ray to either current family cannot enlarge its mathematical cone. Its correlations and compatible same-block cross terms are already admitted.',
        'cases': results})
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    run()
