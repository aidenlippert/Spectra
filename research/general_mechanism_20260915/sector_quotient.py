"""Exact operator-action quotient on an even-N singlet sector.

This is a physical-sector identity construction, not a ground-energy solver or
an equivalence certificate for a pre-existing, truncated SOS ideal. Inputs are
real rational operators of degree at most three, each with one spatial charge.
No determinants, Hamiltonian coefficients, or optimization maps are required.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache

from experiments.marginal_symbolic import adj, canonical, product
from research.interacting_scaling_20260915.singlet_trace import (
    singlet_dimension, singlet_trace, spatial_charge,
)


def psd_quotient(matrix):
    """Exact LDL and a complete kernel, rejecting non-PSD or inexact input."""
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError('A square matrix is required')
    if any(type(x) not in (int, F) for row in matrix for x in row):
        raise ValueError('Only exact real rational entries are accepted')
    gram = [[F(x) for x in row] for row in matrix]
    if any(gram[i][j] != gram[j][i] for i in range(size) for j in range(i)):
        raise ValueError('The matrix must be exactly symmetric')
    remainder = [row[:] for row in gram]
    pivots, diagonal, columns = [], [], []
    for k in range(size):
        value = remainder[k][k]
        if value < 0:
            raise ValueError('Negative exact LDL pivot')
        if not value:
            if any(remainder[k]) or any(row[k] for row in remainder):
                raise ValueError('A zero PSD diagonal must have a zero row')
            continue
        column = [remainder[i][k]/value for i in range(size)]
        pivots.append(k)
        diagonal.append(value)
        columns.append(column)
        for i in range(k, size):
            for j in range(k, size):
                remainder[i][j] -= value*column[i]*column[j]
    if any(x for row in remainder for x in row):
        raise ValueError('Exact LDL reconstruction failed')
    free = [i for i in range(size) if i not in pivots]
    kernel = []
    coordinates = [[F(int(p == j)) for j in range(size)] for p in pivots]
    for j in free:
        vector = [F(0)]*size
        vector[j] = F(1)
        for p, column in reversed(list(zip(pivots, columns))):
            vector[p] = -sum(column[i]*vector[i] for i in range(p+1, size))
        if any(sum(a*b for a, b in zip(row, vector)) for row in gram):
            raise ValueError('Exact kernel reconstruction failed')
        kernel.append(vector)
        for i, p in enumerate(pivots):
            coordinates[i][j] = -vector[p]
    # O_j P = sum_a coordinates[a][j] O_pivot[a] P. This check
    # independently binds those coordinates to the supplied full Gram matrix.
    sparse_coordinates = [[(a, row[j]) for a, row in enumerate(coordinates) if row[j]]
                          for j in range(size)]
    for i in range(size):
        for j in range(size):
            value = sum(x*gram[pivots[a]][pivots[b]]*y
                        for a, x in sparse_coordinates[i]
                        for b, y in sparse_coordinates[j])
            if value != gram[i][j]:
                raise ValueError('Quotient Gram reconstruction failed')
    return {'pivots': pivots, 'positive_diagonal': diagonal,
            'ldl_columns': columns, 'kernel': kernel, 'coordinates': coordinates}


def operator_quotient(operators, modes, particles):
    """Find all singlet-annihilating combinations in the supplied dictionary.

    Different spatial-charge groups are exactly trace-orthogonal. Within each
    group the faithful normalized singlet trace gives a small rational Gram
    matrix. Its kernel is exactly the space of operators annihilating every
    state in the specified singlet sector, not just a trial state.
    """
    dimension = singlet_dimension(modes, particles)
    if not operators:
        raise ValueError('A nonempty operator dictionary is required')
    polys, groups = [], defaultdict(list)
    for index, operator in enumerate(operators):
        for word, value in operator.items():
            if type(value) not in (int, F):
                raise ValueError('Only exact real rational coefficients are accepted')
            if len(word) > 3 or any(type(c) is not int or c not in (0, 1)
                                   or type(i) is not int or not 0 <= i < modes
                                   for c, i in word):
                raise ValueError('Supported input: valid CAR words through degree three')
        poly = canonical(operator)
        charge = spatial_charge(poly) if poly else ()
        groups[charge].append(index)
        polys.append(poly)

    @lru_cache(maxsize=None)
    def trace_word(word):
        return singlet_trace({word: F(1)}, modes, particles)

    records = []
    evaluations = 0
    for charge, indices in sorted(groups.items()):
        size = len(indices)
        gram = [[F(0) for _ in indices] for _ in indices]
        for i, left in enumerate(indices):
            for j in range(i, size):
                terms = product(adj(polys[left]), polys[indices[j]])
                value = sum(c*trace_word(w) for w, c in terms.items())
                gram[i][j] = gram[j][i] = value
                evaluations += 1
        result = psd_quotient(gram)
        records.append({'spatial_charge': charge, 'indices': indices,
                        'gram': gram, **result})
    rank = sum(len(group['pivots']) for group in records)
    return {'modes': modes, 'particles': particles, 'spin': 0,
            'sector_dimension': dimension, 'dictionary_size': len(polys),
            'rank': rank, 'nullity': len(polys)-rank,
            'trace_product_evaluations': evaluations,
            'ungrouped_upper_triangle_pairs': len(polys)*(len(polys)+1)//2,
            'largest_trace_block': max(len(x['indices']) for x in records),
            'distinct_trace_words': trace_word.cache_info().currsize,
            'determinants_enumerated': 0, 'groups': records,
            'existing_SOS_ideal_equivalence_proved': False,
            'ground_energy_improvement_proved': False}
