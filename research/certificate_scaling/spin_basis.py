"""Exact highest-weight decomposition of bounded-degree CAR word spans.

Local kernels use fixed spatial creation/annihilation counts, never occupation
states. Complete descendant rank and spin closure are checked for every pattern.
"""
from fractions import Fraction as F
from math import comb,lcm
from collections import defaultdict
from experiments.marginal_symbolic import canonical,mono,add,scale,validate_word


def spin_weight2(w):
    return sum((2*c-1)*(1 if i%2==0 else -1) for c,i in w)


def spin_action(p,raising=True):
    out={}
    for w,v in p.items():
        for k,(c,i) in enumerate(w):
            active=(c==1 and i%2==(1 if raising else 0)) or (c==0 and i%2==(0 if raising else 1))
            if not active:continue
            replacement=(c,i^1);word=w[:k]+(replacement,)+w[k+1:]
            out=add(out,scale(canonical(mono(word)),v*(1 if c else -1)))
    return out


def rref(rows,ncols):
    A=[[F(x) for x in row] for row in rows];pivots=[];r=0
    for c in range(ncols):
        pivot=next((i for i in range(r,len(A)) if A[i][c]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r];q=A[r][c];A[r]=[x/q for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][c]:
                q=A[i][c];A[i]=[x-q*y for x,y in zip(A[i],A[r])]
        pivots.append(c);r+=1
        if r==len(A):break
    return A,pivots


def nullspace(rows,ncols):
    A,pivots=rref(rows,ncols);vectors=[]
    for free in range(ncols):
        if free in pivots:continue
        v=[F(0)]*ncols;v[free]=F(1)
        for row,p in zip(A,pivots):v[p]=-row[free]
        vectors.append(v)
    return vectors


def decompose_words(words,modes,parity_masks=()):
    import time
    start=time.monotonic()
    if type(modes) is not int or modes<2 or modes%2:raise ValueError('Paired spin modes required')
    for mask in parity_masks:
        if type(mask) is not int or not 0<=mask<1<<modes or any(((mask>>(2*i))&1)!=((mask>>(2*i+1))&1) for i in range(modes//2)):
            raise ValueError('Parity must commute with spin')
    base=set()
    for w in words:
        w=validate_word(w,modes,3);p=canonical(mono(w))
        if len(p)!=1:raise ValueError('Each input word must normalize to one nonzero monomial')
        base.add(next(iter(p)))
    patterns=defaultdict(list)
    for w in sorted(base):
        key=(tuple(sorted(i//2 for c,i in w if c)),tuple(sorted(i//2 for c,i in w if not c)))
        patterns[key].append(w)
    gathered=defaultdict(list);max_dimension=0;highest_count=0;descendant_count=0;actions=0
    def act(p,raising=True):
        nonlocal actions
        actions+=sum(map(len,p))
        return spin_action(p,raising)
    for key,local in sorted(patterns.items()):
        localset=set(local);max_dimension=max(max_dimension,len(local));descendants=[]
        for w in local:
            for raising in (True,False):
                action=act({w:F(1)},raising)
                if not set(action)<=localset:raise ValueError('Input span is not closed under spin')
        for j2 in range(4):
            candidates=[w for w in local if spin_weight2(w)==j2]
            if not candidates:continue
            images=[act({w:F(1)}) for w in candidates]
            imagewords=sorted({w for p in images for w in p})
            rows=[[p.get(w,F(0)) for p in images] for w in imagewords]
            for v in nullspace(rows,len(candidates)):
                highest={w:c for w,c in zip(candidates,v) if c};ladder=[highest]
                for r in range(1,j2+1):
                    ladder.append(scale(act(ladder[-1],False),F(1,r)))
                if act(ladder[-1],False):raise AssertionError('Bad terminal weight')
                for r,p in enumerate(ladder):
                    if not p or not set(p)<=localset or any(spin_weight2(w)!=j2-2*r for w in p):raise AssertionError('Bad descendant')
                    if r and act(p)!=scale(ladder[r-1],j2-r+1):raise AssertionError('Bad raising coefficient')
                word=next(iter(highest));charge=sum(2*c-1 for c,_ in word)
                signature=tuple(sum((mask>>i)&1 for c,i in word)%2 for mask in parity_masks)
                gathered[(charge,j2,signature)].append(ladder)
                highest_count+=1;descendant_count+=len(ladder);descendants.extend(ladder)
        # Exact spanning/injectivity gate in a constant-size local space.
        matrix=[[p.get(w,F(0)) for p in descendants] for w in local]
        _,pivots=rref(matrix,len(descendants))
        if len(descendants)!=len(local) or len(pivots)!=len(local):raise AssertionError('Incomplete or dependent spin decomposition')
    groups=[]
    for (charge,j2,signature),copies in sorted(gathered.items()):
        norm=lcm(*(comb(j2,r) for r in range(j2+1)))
        groups.append({'charge':charge,'two_spin':j2,'signature':signature,'copies':copies,
                       'weights':[norm//comb(j2,r) for r in range(j2+1)]})
    return groups,{'original_dimension':len(base),'highest_weight_dimension':highest_count,
                   'descendant_dimension':descendant_count,'local_pattern_count':len(patterns),
                   'max_local_dimension':max_dimension,'ladder_word_letter_visits':actions,
                   'invariant_Gram_entries':sum(len(g['copies'])**2 for g in groups),
                   'seconds':time.monotonic()-start,'scope':'Exact word-span decomposition and invariant Gram dimensions; no accuracy or solver-speed guarantee.'}
