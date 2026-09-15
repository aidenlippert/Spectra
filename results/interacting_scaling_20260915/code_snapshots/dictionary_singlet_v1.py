"""Direct spin-closed local and collective dictionaries; no global cubic pairs."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from experiments.marginal_symbolic import add, adj, canonical, mono, multiplier_basis, scale


def support(word):
    return frozenset(i//2 for c, i in word)


def weight(word):
    return sum((2*c-1)*(1 if i % 2 == 0 else -1) for c, i in word)


def windows(spatial, width):
    if type(width) is not int or not 2 <= width <= spatial:
        raise ValueError('Cluster width must be between two and the orbital count')
    return [tuple(range(start, start+width)) for start in range(spatial-width+1)]


def clustered_frame(modes, clusters, collective_pairs=True, complete=False):
    import numpy as np
    from research.transfer_solver_20260915.algebra import highest_weights
    if type(modes) is not int or modes < 4 or modes % 2:
        raise ValueError('Paired spin orbitals required')
    spatial = modes//2
    clusters = [tuple(c) for c in clusters]
    if any(len(c) < 2 or len(set(c)) != len(c) or any(type(i) is not int or not 0 <= i < spatial for i in c) for c in clusters):
        raise ValueError('Invalid declared spatial cluster')
    groups, specs = [], []

    def append_block(name, families, spin_reduce):
        local = []
        for family, words in families:
            parts = defaultdict(list)
            for word in dict.fromkeys(words):
                charge = sum(2*c-1 for c, i in word)
                parts[(charge, weight(word))].append(word)
            for signature, entries in sorted(parts.items()):
                local.append({'name': name+':'+family+str(signature),
                    'family': family, 'signature': signature, 'words': entries})
        if spin_reduce:
            local_specs = highest_weights(local)
        else:
            local_specs = [(i, np.eye(len(g['words'])), {'kind': 'full_quadratic'}) for i, g in enumerate(local)]
        offset = len(groups)
        groups.extend(local)
        specs.extend((offset+i, V, {**details, 'declared_block': name}) for i, V, details in local_specs)

    def dagger(w):
        return tuple((1-c, i) for c, i in reversed(w))

    linear = [((0, i),) for i in range(modes)]
    pairs = [((0, j), (0, i)) for i, j in combinations(range(modes), 2)]
    ph = [((1, i), (0, j)) for i in range(modes) for j in range(modes)]
    append_block('global_quadratic', [('linear-', linear), ('linear+', list(map(dagger, linear))),
        ('pair-', pairs), ('pair+', list(map(dagger, pairs))), ('particle-hole', ph)], False)

    def cubic(name, spatial_supports):
        words, triples, linear_words = set(), set(), set()
        for positions in spatial_supports:
            orbitals = tuple(2*p+spin for p in positions for spin in (0, 1))
            linear_words.update(((0, i),) for i in orbitals)
            for i, j in combinations(orbitals, 2):
                words.update(((1, k), (0, j), (0, i)) for k in orbitals)
            triples.update(tuple((0, i) for i in triple) for triple in combinations(orbitals, 3))
        mixed = sorted(linear_words | words, key=lambda w: (len(w), w))
        triple_words = sorted(triples)
        append_block(name, [('mixed-', mixed), ('mixed+', list(map(dagger, mixed))),
            ('triples-', triple_words), ('triples+', list(map(dagger, triple_words)))], True)

    if complete:
        cubic('complete_cubic', [tuple(range(spatial))])
    else:
        if collective_pairs:
            # All pair-supported words share a Gram block per spin channel.
            # Cross terms between distant orbital pairs are retained, unlike
            # a sum of independently optimized two-orbital fragment energies.
            cubic('collective_pair_supports', combinations(range(spatial), 2))
        for index, cluster in enumerate(clusters):
            cubic(f'cluster_{index}', [cluster])
    return groups, specs


def ideal_basis(modes, clusters, complete=False):
    """Global one-body freedoms and declared two-body freedoms, with no local N."""
    polynomials = {}

    def keep(poly):
        poly = canonical(poly)
        if any(weight(w) for w in poly):
            return
        key = tuple(sorted(poly.items()))
        polynomials[key] = poly

    for p in multiplier_basis(modes, max_body=1):
        keep(p)
    if complete:
        for p in multiplier_basis(modes, max_body=2):
            keep(p)
    else:
        # Long-range occupation products are collective number freedoms.
        for i, j in combinations(range(modes), 2):
            keep(mono(((1, i), (1, j), (0, j), (0, i))))
        for cluster in clusters:
            mapping = [2*p+spin for p in cluster for spin in (0, 1)]
            for p in multiplier_basis(len(mapping), max_body=2):
                keep({tuple((c, mapping[i]) for c, i in w): v for w, v in p.items()})
    return list(polynomials.values())


def representative(word):
    """One member of a real Hermitian number-conserving coefficient pair."""
    creators = tuple(i for c, i in word if c)
    annihilators = tuple(i for c, i in word if not c)
    if len(creators) != len(annihilators):
        raise ValueError('A coefficient row must conserve global particle number')
    if creators > annihilators:
        creators, annihilators = annihilators, creators
    return tuple((1, i) for i in creators)+tuple((0, i) for i in annihilators)


def close_rows(words):
    from research.sector_quotient_20260914.fast_twirl import twirl
    rows = {representative(w) for w in words}
    pending = list(rows)
    while pending:
        w = pending.pop()
        p = canonical(mono(w))
        q = twirl(add(p, canonical(adj(p))))
        for v in q:
            r = representative(v)
            if r not in rows:
                rows.add(r)
                pending.append(r)
    return sorted(rows, key=lambda w: (len(w), w))
