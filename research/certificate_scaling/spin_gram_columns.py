"""Exact spin Gram columns with cached canonical adjoints and fused accumulation."""
from functools import lru_cache
from experiments.marginal_symbolic import word_product


@lru_cache(maxsize=100000)
def adjoint_word(w):
    creators=tuple(i for c,i in w if c==1)
    annihilators=tuple(i for c,i in w if c==0)
    expected=tuple((1,i) for i in sorted(set(creators)))+tuple((0,i) for i in sorted(set(annihilators)))
    if w!=expected:raise ValueError('Canonical CAR monomial required')
    parity=(len(creators)*(len(creators)-1)+len(annihilators)*(len(annihilators)-1))//2
    return tuple((1,i) for i in annihilators)+tuple((0,i) for i in creators),(-1 if parity%2 else 1)


def adjoint(p):
    out={}
    for w,c in p.items():
        v,s=adjoint_word(w);out[v]=c if s==1 else -c
    return out


def columns(group):
    copies=group['copies'];weights=group['weights'];cols=[];ix=[];work=0
    if any(len(copy)!=len(weights) for copy in copies):raise ValueError('Channel count mismatch')
    lefts=[[[ (w,c*weight) for w,c in adjoint(p).items()] for p,weight in zip(copy,weights)] for copy in copies]
    for i,left in enumerate(copies):
        for j in range(i,len(copies)):
            result={}
            for channel,(p,q) in enumerate(zip(left,copies[j])):
                for u,a in lefts[i][channel]:
                    for v,b in q.items():
                        ab=a*b
                        for w,s in word_product(u,v):
                            result[w]=result.get(w,0)+(ab if s==1 else -ab)
                work+=len(p)*len(q)
            result={w:c for w,c in result.items() if c}
            if i!=j:
                for w,c in adjoint(result).items():result[w]=result.get(w,0)+c
                result={w:c for w,c in result.items() if c}
            if result!=adjoint(result):raise AssertionError('Non-Hermitian invariant column')
            cols.append(result);ix.append(i*len(copies)+j)
    return cols,ix,work
