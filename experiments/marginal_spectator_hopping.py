"""Fourteen PH-even, reflection-odd one-spectator hopping telescopes.

C(i,j,k)=q_k^power B(i,j), with distinct spectator k. Odd-distance
hopping uses power2; even-distance hopping uses power1. Y=C-C_reflected
on five sites, and T=Y_left-Y_right on six sites cancels under translation.
"""
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_transfer_verify import apply_word


LABELS = {f'{i},{j},{k}': (i,j,k,2 if (j-i)%2 else 1)
          for i in range(5) for j in range(i+1,5) for k in range(5)
          if k not in (i,j) and (i,j,k) < (4-j,4-i,4-k)}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 14:
        raise ValueError('One through fourteen spectator hopping components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical distinct spectator hopping triple required')
    result = {key: _exact(value) for key,value in source.items()}
    result = {key: value for key,value in result.items() if value}
    if not result:
        raise ValueError('Nonzero exact spectator hopping telescope required')
    return result


def actions(terms):
    pairs = {}
    for label,value in terms.items():
        i,j,k,power = LABELS[label]
        for left,right,spectator,sign in [(i,j,k,1),(4-j,4-i,4-k,-1),
                                          (i+1,j+1,k+1,-1),(5-j,5-i,5-k,1)]:
            weights = pairs.setdefault((left,right),{})
            key = (spectator,power)
            weights[key] = weights.get(key,0)+sign*value
    words = {pair: [word for spin in (0,1)
                   for word in [((1,2*pair[0]+spin),(0,2*pair[1]+spin)),
                                ((1,2*pair[1]+spin),(0,2*pair[0]+spin))]]
             for pair in pairs}
    result = []
    for state in range(4096):
        charges = [((state>>(2*i))&3).bit_count()-1 for i in range(6)]
        image = {}
        for pair,weights in pairs.items():
            value = sum(a*charges[k]**power for (k,power),a in weights.items())
            if not value:
                continue
            for word in words[pair]:
                target = apply_word(word,state)
                if target is not None:
                    s,phase = target
                    image[s] = image.get(s,0)+value*phase
        result.append({s:a for s,a in image.items() if a})
    for state,image in enumerate(result):
        reflected,phase = _reflection(state,6)
        reflected_image = {}
        for target,a in image.items():
            if _sector(target,6) != _sector(state,6):
                raise ValueError('Spectator hopping leaves its spin sector')
            if result[target].get(state,0) != a:
                raise ValueError('Spectator hopping is not Hermitian')
            r,sign = _reflection(target,6)
            reflected_image[r] = sign*a
        if reflected_image != {s:phase*a for s,a in result[reflected].items()}:
            raise ValueError('Spectator hopping does not preserve reflection')
    return result
