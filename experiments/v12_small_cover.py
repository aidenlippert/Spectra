"""Untrusted, bounded overlap proposals with common dyadic split weights.

The unchanged V8 checker establishes every accepted cover. This constructor
uses only already constructed ordinary partitions and one coordinate sweep;
it does not claim convergence, novelty or autonomous acquisition.
"""
from fractions import Fraction as F
from math import sqrt,isfinite
from time import perf_counter
from .v7_certificate import clean,BudgetExceeded
from .v8_fractional_cover import _witness


def small_cover(op,partitions,grid_bits=20):
    start=perf_counter()
    if type(grid_bits)is not int or not 1<=grid_bits<=32:raise ValueError('dyadic grid budget')
    if not isinstance(op,dict) or len(op)>40:raise BudgetExceeded('small-cover term budget')
    if not op:return _witness([]),dict(groups=0,incidences=0,construction_seconds=perf_counter()-start)
    n=len(next(iter(op)));op=clean(op,n,40)
    if not op:return _witness([]),dict(groups=0,incidences=0,construction_seconds=perf_counter()-start)
    if not isinstance(partitions,dict) or set(partitions)!={'weighted','firstfit'}:raise ValueError('two ordinary partitions required')
    family=[];known=set()
    for mode in ('weighted','firstfit'):
        seen=[]
        for group in partitions[mode]:
            labels=tuple(sorted(group['labels']))
            if not labels or len(set(labels))!=len(labels) or any(p not in op for p in labels):raise ValueError('partition labels')
            seen.extend(labels)
            if labels not in known:family.append(labels);known.add(labels)
        if len(seen)!=len(op) or set(seen)!=set(op):raise ValueError('not an exact partition')
    if len(family)>80:raise BudgetExceeded('small-cover group budget')
    incidence={p:[] for p in op}
    for j,group in enumerate(family):
        for p in group:incidence[p].append(j)
    if any(not 1<=len(js)<=2 for js in incidence.values()):raise ValueError('incidence multiplicity')
    scale=max(map(abs,op.values()))
    vals={p:float(c/scale) for p,c in op.items()}
    x=[{p:vals[p]/len(incidence[p]) for p in group} for group in family]
    squares=[sum(v*v for v in row.values()) for row in x]
    updates=roots=0
    for p in sorted(op):
        js=incidence[p]
        if len(js)<2:continue
        # Incrementally maintained squares are untrusted numerical work. The
        # exact exported witness covers rounding/cancellation mistakes here.
        others=[sqrt(max(0.,squares[j]-x[j][p]**2)) for j in js]
        roots+=len(js);total=sum(others)
        if total:
            for j,a in zip(js,others):
                old=x[j][p];new=vals[p]*a/total
                if not isfinite(new):raise ValueError('nonfinite split proposal')
                x[j][p]=new;squares[j]+=new*new-old*old;updates+=1
    M=1<<grid_bits;splits=[{} for _ in family];exports=0
    for p,js in incidence.items():
        scores=[abs(x[j][p]) for j in js]
        if any(not isfinite(s) for s in scores):raise ValueError('nonfinite weights')
        total=sum(scores)
        if not total:scores=[1.]*len(js);total=float(len(js))
        anchor=max(range(len(js)),key=lambda i:scores[i])
        remaining=M
        for i,j in enumerate(js):
            if i==anchor:continue
            a=min(remaining,max(0,int(M*scores[i]/total)))
            remaining-=a
            if a:splits[j][p]=op[p]*F(a,M);exports+=1
        if remaining:splits[js[anchor]][p]=op[p]*F(remaining,M);exports+=1
    witness=_witness(splits)
    return witness,dict(groups=len(family),incidences=sum(map(len,family)),coordinate_updates=updates,
                        numerical_roots=roots,dyadic_exports=exports,grid_bits=grid_bits,
                        construction_seconds=perf_counter()-start)


def evolution_witness(cover,degree,derivative_nonzero,weighted_bound):
    empty=_witness([])
    if not derivative_nonzero:ws={'jump:0':empty,'residual:0:0':empty}
    else:
        ws={'jump:0':empty,**{f'residual:0:{k}':empty for k in range(degree+1)}}
        ws[f'residual:0:{degree}']=dict(cover,groups=[
            dict(coefficients={p:str(-F(c)) for p,c in g['coefficients'].items()},upper=g['upper'])
            for g in cover['groups']])
    return dict(schema='v8-evolution-cover-1',witnesses=ws,claimed_bound=str(weighted_bound))
