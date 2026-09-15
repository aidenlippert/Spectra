"""Trajectory-wide residual-driven coordinate Galerkin; supplied baseline math."""
from fractions import Fraction as F
from time import perf_counter
from .v7_certificate import (Generator,Piece,clean,rational,norm_witness,check_certificate,BudgetExceeded)


def projected_taylor(gen,initial,basis,duration,max_order=24):
    """Fixed-space recurrence with exact omitted derivatives, including tail."""
    if type(max_order) is not int or not 0<=max_order<=24:raise ValueError('order cap')
    if rational(duration)<=0:raise ValueError('duration')
    B=set(clean({p:F(1) for p in basis},gen.n,gen.max_terms))
    c={p:v for p,v in clean(initial,gen.n,gen.max_terms).items() if p in B}
    cs=[c];omitted=[]
    for k in range(max_order+1):
        full=gen.apply(c);omitted.append({p:v for p,v in full.items() if p not in B})
        if k<max_order:
            c={p:v/F(k+1) for p,v in full.items() if p in B};cs.append(c)
    return tuple(cs),tuple(omitted)


def _witness(omitted,tail,duration,tolerance):
    # Lower derivatives cancel only within the fixed Galerkin basis.
    residuals=[{p:-v for p,v in c.items()} for c in omitted[:-1]]+[{p:-v for p,v in tail.items()}]
    while len(residuals)>1 and not residuals[-1]:residuals.pop()
    best=None;cost=dict(group_comparisons=0,sorted_terms=0,sqrt_enclosures=0)
    for mode in ('l1','firstfit','weighted'):
        ws={'jump:0':[]};bound=F(0)
        for k,op in enumerate(residuals):
            groups,counts=norm_witness(op,mode);ws[f'residual:0:{k}']=groups
            bound+=duration**(k+1)/F(k+1)*sum((F(g['upper']) for g in groups),F(0))
            for key in cost:cost[key]+=counts[key]
        w=dict(schema='v7-residual-1',integration_basis='power',witnesses=ws,claimed_bound=str(bound),construction_cost=dict(cost))
        if best is None or bound<F(best['claimed_bound']):best=w
        if bound<=tolerance:break
    return best,cost


def adaptive_galerkin(gen,initial,duration,tolerance,*,max_order=24,max_expansions=16,batch=None,growth='residual'):
    """Increment order on each fixed space; enrich from accumulated residual.

An exact l1 bound can accept early. Otherwise enrich once the in-space tail
is <= eps/4, using accumulated omitted coefficients over the computed trajectory.
These are search choices, not necessity or impossibility certificates.
"""
    start=perf_counter();T=rational(duration);eps=rational(tolerance)
    if T<=0 or eps<0:raise ValueError('duration/tolerance')
    if type(max_order) is not int or not 0<=max_order<=24:raise ValueError('order cap')
    if type(max_expansions) is not int or not 0<=max_expansions<=16:raise ValueError('expansion cap')
    if batch is not None and (type(batch) is not int or not 1<=batch<=512):raise ValueError('batch cap')
    if growth not in ('residual','bfs'):raise ValueError('growth policy')
    seed=clean(initial,gen.n,gen.max_terms);B=set(seed);rounds=[];checks=[];norm_cost=dict(group_comparisons=0,sorted_terms=0,sqrt_enclosures=0)
    frontier=set(seed)
    def receipt(status,reason=None,**extra):
        return dict(status=status,reason=reason,basis=tuple(sorted(B)),rounds=rounds,checker_receipts=checks,
                    generator_cost=dict(gen.cost),norm_work=dict(norm_cost),complete_seconds=perf_counter()-start,**extra)
    try:
        for expansion in range(max_expansions+1):
            cs=[dict(seed)];omitted=[];scores={};omission_l1=F(0)
            for m in range(max_order+1):
                full=gen.apply(cs[-1]);factor=T**(m+1)/F(m+1)
                outside={p:v for p,v in full.items() if p not in B};omitted.append(outside)
                inside_l1=factor*sum((abs(v) for p,v in full.items() if p in B),F(0))
                outside_l1=factor*sum((abs(v) for v in outside.values()),F(0))
                omission_l1+=outside_l1
                for p,v in outside.items():scores[p]=scores.get(p,F(0))+factor*abs(v)
                total_l1=omission_l1+inside_l1
                if total_l1<=eps or inside_l1<=eps/4 or m==max_order:
                    break
                cs.append({p:v/F(m+1) for p,v in full.items() if p in B})
            if sum(map(len,cs))>50000:raise BudgetExceeded('polynomial coefficient budget')
            w,nc=_witness(omitted,full,T,eps)
            for key in norm_cost:norm_cost[key]+=nc[key]
            round_receipt=dict(expansion=expansion,basis_size=len(B),order=m,inside_tail_l1=str(inside_l1),omitted_l1=str(omission_l1),bound=w['claimed_bound'],scored_terms=len(scores))
            rounds.append(round_receipt)
            if F(w['claimed_bound'])<=eps:
                piece=Piece(T,tuple(cs));checker_gen=Generator(gen.h,gen.gamma,gen.n,gen.max_terms)
                t0=perf_counter();ans=check_certificate(checker_gen,seed,[piece],w,eps,expected_time=T)
                checks.append(dict(receipt=ans,checking_seconds=perf_counter()-t0))
                if ans['status']!='certified':return receipt('refused','checker rejected candidate: '+str(ans))
                return receipt('certified',piece=piece,certificate=w,checker=ans,expansions=expansion)
            if expansion==max_expansions:return receipt('refused','expansion budget')
            available=gen.max_terms-len(B)
            if not available:return receipt('refused','basis term budget')
            if growth=='bfs':
                neighbors=set()
                for p in sorted(frontier):neighbors.update(gen.apply({p:F(1)}))
                candidates=sorted(neighbors-B)
            else:
                candidates=sorted((p for p in scores if p not in B),key=lambda p:(-scores[p],p))
            if not candidates:return receipt('refused','order budget without omitted frontier')
            if growth=='bfs' and len(candidates)>available:return receipt('refused','basis term budget')
            amount=len(candidates) if growth=='bfs' else min(batch if batch is not None else max(4,len(B)//2),available)
            new=candidates[:amount];round_receipt['added_labels']=new;B.update(new);frontier=set(new)
        return receipt('refused','expansion budget')
    except (ValueError,RuntimeError) as exc:
        return receipt('refused',str(exc))
