"""Exact range-three hopping telescope, with no physical range-three coupling.

B(i,j) sums both spin hoppings. Y=B(0,3)-B(1,4) on five sites;
T=Y_left-Y_right=B(0,3)-2B(1,4)+B(2,5) on six sites.
"""
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_transfer_verify import apply_word


def coefficient(value):
    value = _exact(value)
    if not value:
        raise ValueError('Nonzero exact hopping telescope required')
    return value


def actions():
    words = []
    for i, weight in enumerate((1, -2, 1)):
        for spin in (0, 1):
            a, b = 2*i+spin, 2*(i+3)+spin
            words.extend([(((1, a), (0, b)), weight), (((1, b), (0, a)), weight)])
    result = []
    for state in range(4096):
        image = {}
        for word, weight in words:
            target = apply_word(word, state)
            if target is not None:
                s, phase = target
                image[s] = image.get(s, 0)+weight*phase
        result.append({s: a for s, a in image.items() if a})
    for state, image in enumerate(result):
        reflected, phase = _reflection(state, 6)
        reflected_image = {}
        for target, a in image.items():
            if _sector(target, 6) != _sector(state, 6):
                raise ValueError('Hopping telescope leaves its spin sector')
            if result[target].get(state, 0) != a:
                raise ValueError('Hopping telescope is not Hermitian')
            r, sign = _reflection(target, 6)
            reflected_image[r] = sign*a
        if reflected_image != {s: phase*a for s, a in result[reflected].items()}:
            raise ValueError('Hopping telescope does not preserve reflection')
    return result


def projected_matrix(columns, action):
    """E^T T E in the physical reflection basis used by local positivity."""
    lookup = {s: (i, a) for i, col in enumerate(columns) for s, a in col.items()}
    matrix = [[0 for _ in columns] for _ in columns]
    for j, col in enumerate(columns):
        for s, a in col.items():
            for target, b in action[s].items():
                if target in lookup:
                    i, c = lookup[target]
                    matrix[i][j] += a*b*c
    return matrix
