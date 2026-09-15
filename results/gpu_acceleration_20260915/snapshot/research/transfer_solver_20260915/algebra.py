"""Input-derived spin dictionaries for the existing exact SOS checker.

The rational highest-weight kernel is the supplied H8 construction, with group
indices replaced by charge/spin/parity metadata and mixed linear words retained
uniformly. No saved factor, Gram matrix, coefficient map, or MPS is an input.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
import numpy as np
from scipy import sparse
from experiments.marginal_coefficient import dictionaries
from experiments.marginal_molecule_stress import parity_generators
from experiments.marginal_symbolic import add, mono, multiplier_basis


def weight(word):
    return sum((2*c-1)*(1 if i % 2 == 0 else -1) for c, i in word)


def normalize(word):
    cr = [i for c, i in word if c]
    an = [i for c, i in word if not c]
    if len(set(cr)) != len(cr) or len(set(an)) != len(an):
        return None, 0
    sign = (-1)**(sum(cr[i] > cr[j] for i in range(len(cr)) for j in range(i+1, len(cr)))
                  + sum(an[i] < an[j] for i in range(len(an)) for j in range(i+1, len(an))))
    return tuple([(1, i) for i in sorted(cr)]+[(0, i) for i in sorted(an, reverse=True)]), sign


def frame(h, modes):
    # Restrict binary phase symmetries to those commuting with spin rotations.
    # Added ladder words constrain only this support analysis, never physical H.
    support = add(h, *(mono(((1, i), (0, i+1))) for i in range(0, modes, 2)))
    masks = parity_generators(support, modes)
    if any(((mask >> i) & 1) != ((mask >> (i+1)) & 1) for mask in masks for i in range(0, modes, 2)):
        raise AssertionError('Spatial parity must commute with spin raising')
    def signature(word):
        return (sum(2*c-1 for c, i in word), weight(word)) + tuple(
            sum((mask >> i) & 1 for c, i in word) % 2 for mask in masks)
    groups = []
    for block in dictionaries(modes, 'mixed'):
        parts = defaultdict(list)
        for word in block['words']:
            parts[signature(word)].append(word)
        for sig, words in sorted(parts.items()):
            groups.append({'name': block['name']+str(sig), 'family': block['name'],
                           'signature': sig, 'words': words})
    rows = []
    for degree in range(4):
        sets = list(combinations(range(modes), degree))
        for i, a in enumerate(sets):
            for b in sets[i:]:
                word = tuple((1, k) for k in a)+tuple((0, k) for k in b)
                if signature(word) == signature(()):
                    rows.append(word)
    basis = [p for p in multiplier_basis(modes, max_body=2)
             if all(signature(w) == signature(()) for w in p)]
    return groups, rows, basis, {'spatial_parity_masks': masks, 'mixed_linear_words_retained': True,
                                'operator_dictionary_dimensions': [len(g['words']) for g in groups]}


def raising(groups):
    ids = [i for i, g in enumerate(groups) if g['family'] in ('mixed-', 'mixed+')]
    lookup = {tuple(map(tuple, w)): (g, i) for g in ids for i, w in enumerate(groups[g]['words'])}
    # Raising closure alone does not detect deletion of a lower spin-half
    # component. Both ladder directions must act within the declared dictionary.
    for g in ids:
        for word in groups[g]['words']:
            for k, (c, i) in enumerate(word):
                if (c and i % 2 == 0) or (not c and i % 2 == 1):
                    changed = list(word)
                    changed[k] = (c, i ^ 1)
                    w, sign = normalize(changed)
                    if sign and w not in lookup:
                        raise ValueError('Mixed dictionary is not closed under spin lowering')
    matrices, targets = {}, {}
    for g in ids:
        entries = defaultdict(list)
        for j, word in enumerate(groups[g]['words']):
            for k, (c, i) in enumerate(word):
                if (c and i % 2 == 1) or (not c and i % 2 == 0):
                    changed = list(word)
                    changed[k] = (c, i ^ 1)
                    w, sign = normalize(changed)
                    if not sign:
                        continue
                    if w not in lookup:
                        raise ValueError('Mixed dictionary is not closed under spin raising')
                    target, a = lookup[w]
                    entries[target].append((a, j, sign if c else -sign))
        if len(entries) > 1:
            raise ValueError('Spin raising changes the declared spatial symmetry class')
        if entries:
            target, terms = next(iter(entries.items()))
            matrix = sparse.coo_matrix(([v for i, j, v in terms],
                                        ([i for i, j, v in terms], [j for i, j, v in terms])),
                                       shape=(len(groups[target]['words']), len(groups[g]['words'])),
                                       dtype=np.int64).tocsr()
            matrix.sum_duplicates()
            matrix.eliminate_zeros()
            matrices[g], targets[g] = matrix, target
    return ids, matrices, targets


def kernel(U):
    defect = U@U.T-3*sparse.eye(U.shape[0], dtype=np.int64)
    defect.eliminate_zeros()
    if defect.nnz or max(np.asarray(U.getnnz(axis=0)).ravel(), default=0) > 1:
        raise ValueError('Unsupported highest-weight raising-map structure')
    cols = []
    for j in np.flatnonzero(np.asarray(U.getnnz(axis=0)).ravel() == 0):
        v = np.zeros(U.shape[1], dtype=np.int64)
        v[j] = 6
        cols.append(v)
    for i in range(U.shape[0]):
        js = U.indices[U.indptr[i]:U.indptr[i+1]]
        signs = U.data[U.indptr[i]:U.indptr[i+1]]
        if len(js) != 3 or any(abs(s) != 1 for s in signs):
            raise ValueError('Expected three signed kernel entries')
        v = np.zeros(U.shape[1], dtype=np.int64)
        v[js[0]], v[js[1]] = 3*signs[1], -3*signs[0]
        cols.append(v)
        v = np.zeros(U.shape[1], dtype=np.int64)
        v[js[0]], v[js[1]], v[js[2]] = 2*signs[0], 2*signs[1], -4*signs[2]
        cols.append(v)
    K6 = np.array(cols, dtype=np.int64).T
    if np.any(U@K6):
        raise AssertionError('Exact kernel membership failed')
    gram = K6.T@K6
    if np.any(gram-np.diag(np.diag(gram))) or np.any(np.diag(gram) <= 0):
        raise AssertionError('Exact kernel rank failed')
    if K6.shape[1] != U.shape[1]-U.shape[0]:
        raise AssertionError('Incomplete highest-weight kernel')
    return K6/6, {'integer_denominator': 6, 'kernel_columns': K6.shape[1],
                  'integer_kernel_and_rank_verified': True}


def highest_weights(groups):
    ids, matrices, targets = raising(groups)
    result = [(i, np.eye(len(g['words'])), {'kind': 'unchanged_full_dictionary'})
              for i, g in enumerate(groups) if i not in ids]
    covered = set()
    for high in [i for i in ids if i not in targets]:
        chain = [high]
        while any(t == chain[0] for t in targets.values()):
            predecessors = [g for g, t in targets.items() if t == chain[0]]
            if len(predecessors) != 1:
                raise ValueError('Ambiguous spin chain')
            chain.insert(0, predecessors[0])
        if len(chain) != 4 or [weight(groups[g]['words'][0]) for g in chain] != [-3, -1, 1, 3]:
            raise ValueError(('Unsupported mixed spin chain', chain))
        if covered.intersection(chain):
            raise ValueError('Repeated magnetic group')
        covered.update(chain)
        K, receipt = kernel(matrices[chain[2]])
        result.append((chain[2], K, {'kind': 'spin_half_highest', 'chain': chain, **receipt}))
        result.append((high, np.eye(len(groups[high]['words'])),
                       {'kind': 'spin_three_half_highest', 'chain': chain}))
    if covered != set(ids):
        raise ValueError('Incomplete mixed spin decomposition')
    return result
