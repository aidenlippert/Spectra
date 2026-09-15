"""Two purely offdiagonal stationary telescopes after full spin-word closure."""
from collections import defaultdict
from fractions import Fraction as F
from experiments.marginal_local_hubbard_block import _exact
from experiments.marginal_coherent_projector_telescope import build as projector_build, embed

LABELS = {'346,409,1': {346: 1, 409: 1}, '314,614,1': {314: 1, 614: 1}}


def build(label):
    if type(label) is not str or label not in LABELS:
        raise ValueError('Canonical pure-coherence label required')
    denominator, full_y, full_t = projector_build(LABELS[label])
    y = {key: value for key, value in full_y.items() if key[0] != key[1]}
    telescope = {key: value for key, value in full_t.items() if key[0] != key[1]}
    left, right = embed(y, True), embed(y, False)
    expected = {key: left.get(key, 0)-right.get(key, 0) for key in left.keys() | right.keys()}
    if not telescope or telescope != {key: value for key, value in expected.items() if value}:
        raise ValueError('Pure-coherence embedding identity failed')
    # Removing diagonal commutes with every signed permutation checked by projector_build:
    # all inherited Hermiticity, spin, symmetry and periodic identities therefore survive.
    rows = defaultdict(int)
    for (r, c), value in y.items():
        rows[r] += abs(value)
    if F(max(rows.values()), denominator) > F(1, 4):
        raise ValueError('Fixed pure-coherence norm bound failed')
    return denominator, y, telescope


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 2:
        raise ValueError('One or two pure-coherence components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical pure-coherence label required')
    terms = {key: _exact(value) for key, value in source.items()}
    terms = {key: value for key, value in terms.items() if value}
    if not terms:
        raise ValueError('Nonzero exact pure-coherence telescope required')
    return terms


def actions(terms):
    result = [{} for _ in range(4096)]
    for label, coefficient in terms.items():
        denominator, _, matrix = build(label)
        for (row, column), value in matrix.items():
            image = result[column]
            image[row] = image.get(row, F(0)) + coefficient*F(value, denominator)
    return [{s: a for s, a in image.items() if a} for image in result]
