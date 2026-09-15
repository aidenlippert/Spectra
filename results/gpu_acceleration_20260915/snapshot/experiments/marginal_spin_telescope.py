"""Four compact reflection-odd spin-dot telescopes on five sites.

A(i,j)=4 S_i dot S_j; Y=A(i,j)-A(4-j,4-i); T=Y_left-Y_right.
These auxiliary terms cancel under translation and add no physical coupling.
"""
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_transfer_verify import apply_word


LABELS = {'0,1': (0,1), '0,2': (0,2), '0,3': (0,3), '1,2': (1,2)}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 4:
        raise ValueError('One through four spin telescope components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical reflection-odd spin pair required')
    result = {key: _exact(value) for key, value in source.items()}
    result = {key: value for key, value in result.items() if value}
    if not result:
        raise ValueError('Nonzero exact spin telescope required')
    return result


def actions(terms):
    pairs = {}
    for label, value in terms.items():
        i, j = LABELS[label]
        for pair, sign in [((i,j),1), ((4-j,4-i),-1),
                           ((i+1,j+1),-1), ((5-j,5-i),1)]:
            pairs[pair] = pairs.get(pair, 0)+sign*value
    words = []
    for (i,j), value in pairs.items():
        if not value:
            continue
        for a in (0,1):
            for b in (0,1):
                m, n = 2*i+a, 2*j+b
                words.append((((1,m),(0,m),(1,n),(0,n)), value*(-1)**(a+b)))
            words.append((((1,2*i+a),(0,2*i+1-a),
                           (1,2*j+1-a),(0,2*j+a)), 2*value))
    result = []
    for state in range(4096):
        image = {}
        for word, value in words:
            target = apply_word(word, state)
            if target is not None:
                s, phase = target
                image[s] = image.get(s, 0)+value*phase
        result.append({s: a for s, a in image.items() if a})
    for state, image in enumerate(result):
        reflected, phase = _reflection(state, 6)
        reflected_image = {}
        for target, a in image.items():
            if _sector(target, 6) != _sector(state, 6):
                raise ValueError('Spin telescope leaves its spin sector')
            if result[target].get(state, 0) != a:
                raise ValueError('Spin telescope is not Hermitian')
            r, sign = _reflection(target, 6)
            reflected_image[r] = sign*a
        if reflected_image != {s: phase*a for s, a in result[reflected].items()}:
            raise ValueError('Spin telescope does not preserve reflection')
    return result
