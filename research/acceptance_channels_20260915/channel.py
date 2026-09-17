"""Reconstruct the user-supplied H8 channel; no source dual is inferred."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
from experiments.marginal_symbolic import add, adj, canonical, encode, mono, product, scale


def supplied_channel():
    # Even mode labels are alpha; odd labels are beta.
    B = add(
        mono(((1, 4), (0, 5), (0, 7)), F(7423389, 10**7)),
        mono(((1, 4), (0, 7), (0, 9)), F(-6260658, 10**7)),
        mono(((1, 2), (0, 3), (0, 7)), F(2055055, 10**7)),
        mono(((1, 6), (0, 3), (0, 7)), F(1214168, 10**7)),
    )
    A = add(
        mono(((1, 4), (0, 5)), F(7423389, 10**7)),
        mono(((1, 4), (0, 9)), F(6260658, 10**7)),
        mono(((1, 2), (0, 3)), F(2055055, 10**7)),
        mono(((1, 6), (0, 3)), F(1214168, 10**7)),
    )
    P = add(product(adj(B), B), product(B, adj(B)))
    return canonical(B), A, P


def check():
    B, A, P = supplied_channel()
    n = mono(((1, 7), (0, 7)))
    aa = product(A, adj(A))
    commutator = add(product(adj(A), A), scale(aa, -1))
    if canonical(B) != product(A, mono(((0, 7),))):
        raise ValueError('The supplied factorization does not reproduce B')
    if P != add(aa, product(n, commutator)):
        raise ValueError('The occupation-conditioned identity does not reproduce P')
    if canonical(adj(P)) != P:
        raise ValueError('The separator is not Hermitian')
    if len(P) != 31 or max(map(len, P)) != 4:
        raise ValueError('The claimed compact expansion does not reproduce')
    support = sorted({mode // 2 for word in B for _, mode in word})
    windows = [set(range(0, 4)), set(range(2, 6)), set(range(4, 8))]
    return {
        'kind': 'exact_reconstruction_of_supplied_channel',
        'B': encode(B), 'A': encode(A), 'P': encode(P),
        'factorization_verified': True,
        'spatial_support': support,
        'contained_in_any_reported_old_window': any(set(support) <= w for w in windows),
        'P_terms': len(P), 'P_max_degree': max(map(len, P)),
        'commutator_max_degree': max(map(len, commutator)),
        'positivity_argument': 'Exactly B†B + BB† on every Fock state',
        'source_dual_supplied': False,
        'reported_negative_dual_value_replayed': False,
        'new_energy_bound': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = check()
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('A', 'B', 'P')}))
