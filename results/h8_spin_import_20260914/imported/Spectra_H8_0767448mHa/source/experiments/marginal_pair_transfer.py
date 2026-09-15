"""Four reflection-odd pair-transfer differences with zero translated sum.

d_i†=c_i_up† c_i_down†, J(i,j)=d_i†d_j+d_j†d_i.
On five sites Y=J(i,j)-J(4-j,4-i); on six sites T=Y_left-Y_right.
"""
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_transfer_verify import apply_word

LABELS = {'0,1': (0,1), '0,2': (0,2), '0,3': (0,3), '1,2': (1,2)}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 4:
        raise ValueError('One through four pair-transfer components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical pair-transfer label required')
    terms = {key: _exact(value) for key,value in source.items()}
    terms = {key: value for key,value in terms.items() if value}
    if not terms:
        raise ValueError('Nonzero exact pair-transfer telescope required')
    return terms


def actions(terms):
    weights = {}
    for label,value in terms.items():
        i,j = LABELS[label]
        for pair,sign in [((i,j),1), ((4-j,4-i),-1),
                          ((i+1,j+1),-1), ((5-j,5-i),1)]:
            weights[pair] = weights.get(pair,0) + sign*value
    result = []
    for state in range(4096):
        image = {}
        for (i,j),value in weights.items():
            for a,b in ((i,j),(j,i)):
                word = ((1,2*a),(1,2*a+1),(0,2*b+1),(0,2*b))
                target = apply_word(word,state)
                if target is not None:
                    s,phase = target
                    image[s] = image.get(s,0) + value*phase
        result.append({s:a for s,a in image.items() if a})
    for state,image in enumerate(result):
        reflected,phase = _reflection(state,6)
        reflected_image = {}
        for target,a in image.items():
            if _sector(target,6) != _sector(state,6) or (target^state).bit_count() != 4:
                raise ValueError('Pair transfer leaves its spin sector or pair support')
            if result[target].get(state,0) != a:
                raise ValueError('Pair transfer is not Hermitian')
            r,sign = _reflection(target,6)
            reflected_image[r] = sign*a
        if reflected_image != {s:phase*a for s,a in result[reflected].items()}:
            raise ValueError('Pair transfer does not preserve reflection')
    return result
