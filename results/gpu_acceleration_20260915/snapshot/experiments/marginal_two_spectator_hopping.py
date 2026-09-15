"""Thirty PH-even reflected hopping telescopes with two charge spectators."""
from itertools import combinations
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_transfer_verify import apply_word


LABELS = {','.join(map(str,label)):label
          for i in range(5) for j in range(i+1,5)
          for k,l in combinations([s for s in range(5) if s not in (i,j)],2)
          for p,q in (((1,1),(2,2)) if (j-i)%2 else ((1,2),(2,1)))
          for label in [(i,j,k,l,p,q)] if label < (4-j,4-i,4-l,4-k,q,p)}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 30:
        raise ValueError('One through thirty two-spectator hopping components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical two-spectator hopping label required')
    terms = {key:_exact(value) for key,value in source.items()}
    terms = {key:value for key,value in terms.items() if value}
    if not terms:
        raise ValueError('Nonzero exact two-spectator hopping telescope required')
    return terms


def actions(terms):
    pairs = {}
    for label,value in terms.items():
        i,j,k,l,p,q = LABELS[label]
        for left,right,a,b,r,s,sign in [(i,j,k,l,p,q,1),(4-j,4-i,4-l,4-k,q,p,-1),
                                       (i+1,j+1,k+1,l+1,p,q,-1),(5-j,5-i,5-l,5-k,q,p,1)]:
            weights = pairs.setdefault((left,right),{})
            key = (a,b,r,s);weights[key] = weights.get(key,0)+sign*value
    result=[]
    for state in range(4096):
        charges=[((state>>(2*i))&3).bit_count()-1 for i in range(6)]
        image={}
        for (i,j),weights in pairs.items():
            value=sum(a*charges[k]**p*charges[l]**q for (k,l,p,q),a in weights.items())
            if not value:continue
            for left,right in ((i,j),(j,i)):
                for spin in (0,1):
                    target=apply_word(((1,2*left+spin),(0,2*right+spin)),state)
                    if target is not None:
                        s,phase=target;image[s]=image.get(s,0)+value*phase
        result.append({s:a for s,a in image.items() if a})
    for state,image in enumerate(result):
        reflected,phase=_reflection(state,6);reflected_image={}
        for target,a in image.items():
            if _sector(target,6)!=_sector(state,6) or (target^state).bit_count()!=2:
                raise ValueError('Two-spectator hopping leaves its spin sector or hopping support')
            if result[target].get(state,0)!=a:raise ValueError('Two-spectator hopping is not Hermitian')
            r,sign=_reflection(target,6);reflected_image[r]=sign*a
        if reflected_image!={s:phase*a for s,a in result[reflected].items()}:
            raise ValueError('Two-spectator hopping does not preserve reflection')
    return result
