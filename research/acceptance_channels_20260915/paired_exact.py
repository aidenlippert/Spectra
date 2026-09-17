"""Exact quartic evaluation of strictly matched spin-transfer/adjoint factors.

All unmatched blocks use the original expansion. No numerical library is imported.
"""
from collections import defaultdict
from fractions import Fraction
import time
from experiments.marginal_symbolic import expand_squares as original_expand,validate_word,word_product

CALLS=[]


def validated(block,denominator,modes,max_degree):
    if type(denominator) is not int or denominator<=0:raise ValueError('Invalid factor denominator')
    words=[validate_word(w,modes,max_degree) for w in block['words']]
    if not words or len({sum(2*c-1 for c,_ in w) for w in words})!=1:
        raise ValueError('Each factor dictionary must have one charge')
    factors=block['factor']
    if any(len(row)!=len(words) or any(type(c) is not int for c in row) for row in factors):
        raise ValueError('Invalid integer factor')
    return words,factors


def supported(words,modes):
    return modes%2==0 and all(len(w)==3 and [c for c,_ in w]==[1,0,0]
        and w[0][1]%2==0 and w[1][1]%2==1 and w[2][1]%2==1
        and w[1][1]!=w[2][1] for w in words)


def pair_numerator(words,factors):
    total={}
    def emit(word,value):
        if not value:return
        for normal,sign in word_product((),word):
            total[normal]=total.get(normal,0)+sign*value
    for row in factors:
        C={}
        for w,value in zip(words,row):
            p,q,r=[i for _,i in w]
            if q>r:q,r,value=r,q,-value
            key=p,q,r;C[key]=C.get(key,0)+value
        by_pair=defaultdict(list);by_creator=defaultdict(list);by_contraction=defaultdict(list)
        for (p,q,r),value in C.items():
            if not value:continue
            by_pair[q,r].append((p,value));by_creator[p].append((q,r,value))
            by_contraction[r].append((p,q,value));by_contraction[q].append((p,r,-value))
        for items in by_pair.values():
            for p,a in items:
                for s,b in items:emit(((1,p),(0,s)),a*b)
        for items in by_creator.values():
            for q,r,a in items:
                for t,u,b in items:emit(((1,t),(1,u),(0,q),(0,r)),-a*b)
        for items in by_contraction.values():
            for p,w,a in items:
                for s,v,b in items:emit(((1,p),(1,v),(0,w),(0,s)),-a*b)
    return {w:c for w,c in total.items() if c}


def expand_squares(blocks,denominator,modes,max_degree=3):
    started=time.monotonic()
    if type(denominator) is not int or denominator<=0:raise ValueError('Invalid factor denominator')
    checked=[validated(b,denominator,modes,max_degree) for b in blocks]
    fallback=[];integer={};pairs=rows=nonzeros=0;i=0;contraction_seconds=0.
    while i<len(blocks):
        words,factors=checked[i]
        mate=[tuple((1-c,k) for c,k in reversed(w)) for w in words]
        if supported(words,modes) and i+1<len(blocks) and checked[i+1][0]==mate and checked[i+1][1]==factors:
            begin=time.monotonic();addition=pair_numerator(words,factors)
            for w,c in addition.items():integer[w]=integer.get(w,0)+c
            contraction_seconds+=time.monotonic()-begin
            rows+=2*len(factors);nonzeros+=2*sum(c!=0 for row in factors for c in row)
            pairs+=1;i+=2
        else:fallback.append(blocks[i]);i+=1
    begin=time.monotonic()
    polynomial,stats=original_expand(fallback,denominator,modes,max_degree)
    fallback_seconds=time.monotonic()-begin
    for w,c in integer.items():polynomial[w]=polynomial.get(w,Fraction(0))+Fraction(c,denominator**2)
    polynomial={w:c for w,c in polynomial.items() if c}
    stats={'factor_rows':stats['factor_rows']+rows,'factor_nonzeros':stats['factor_nonzeros']+nonzeros}
    CALLS.append({'exact_matched_block_pairs':pairs,'contracted_factor_rows':rows,
        'quartic_contraction_seconds':contraction_seconds,'original_fallback_seconds':fallback_seconds,
        'total_seconds':time.monotonic()-started,'unmatched_blocks':len(fallback)})
    return polynomial,stats
