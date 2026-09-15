"""Eighteen PH-even hopping telescopes with three distinct charge spectators."""
from itertools import combinations, product
from experiments.marginal_local_hubbard_block import _exact, _reflection, _sector
from experiments.marginal_transfer_verify import apply_word


LABELS = {','.join(map(str,label)):label for i,j in combinations(range(5),2)
          for k,l,m in [sorted(set(range(5))-{i,j})]
          for p,q,r in product((1,2),repeat=3) if (p+q+r+j-i)%2==1
          for label in [(i,j,k,l,m,p,q,r)]
          if label < (4-j,4-i,4-m,4-l,4-k,r,q,p)}


def coefficients(source):
    if type(source) is not dict or not 1 <= len(source) <= 18:
        raise ValueError('One through eighteen three-spectator hopping components required')
    if any(type(key) is not str or key not in LABELS for key in source):
        raise ValueError('Canonical three-spectator hopping label required')
    terms = {key:_exact(value) for key,value in source.items()}
    terms = {key:value for key,value in terms.items() if value}
    if not terms:
        raise ValueError('Nonzero exact three-spectator hopping telescope required')
    return terms


def actions(terms):
    pairs = {}
    for label,value in terms.items():
        i,j,k,l,m,p,q,r = LABELS[label]
        mirror=(4-j,4-i,4-m,4-l,4-k,r,q,p)
        shift=lambda t:tuple(s+1 for s in t[:5])+t[5:]
        for left,right,a,b,c,u,v,w,sign in [(*LABELS[label],1),(*mirror,-1),
                                          (*shift(LABELS[label]),-1),(*shift(mirror),1)]:
            weights=pairs.setdefault((left,right),{})
            key=(a,b,c,u,v,w);weights[key]=weights.get(key,0)+sign*value
    result=[]
    for state in range(4096):
        charges=[((state>>(2*i))&3).bit_count()-1 for i in range(6)]
        image={}
        for (i,j),weights in pairs.items():
            value=sum(a*charges[k]**p*charges[l]**q*charges[m]**r
                      for (k,l,m,p,q,r),a in weights.items())
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
                raise ValueError('Three-spectator hopping leaves its spin sector or hopping support')
            if result[target].get(state,0)!=a:raise ValueError('Three-spectator hopping is not Hermitian')
            r,sign=_reflection(target,6);reflected_image[r]=sign*a
        if reflected_image!={s:phase*a for s,a in result[reflected].items()}:
            raise ValueError('Three-spectator hopping does not preserve reflection')
    return result
