"""Bounded overlapping anticommuting coefficient-split certificates (v8)."""
from fractions import Fraction
from math import sqrt
from time import perf_counter
from itertools import combinations

MAX_TERMS=512; MAX_INCIDENCES=2048; MAX_GROUPS=128; MAX_BITS=8192

def anti(a,b):
    if len(a)!=len(b): return False
    odd=0
    for x,y in zip(a,b):
        if x!='I' and y!='I' and x!=y: odd ^= 1
    return bool(odd)

def rat(x):
    if isinstance(x,bool): raise ValueError('boolean rational')
    if isinstance(x,Fraction): q=x
    elif isinstance(x,int): q=Fraction(x)
    elif isinstance(x,str): q=Fraction(x)
    else: raise ValueError('exact rational required')
    if max(q.numerator.bit_length(),q.denominator.bit_length())>MAX_BITS: raise ValueError('rational bit budget')
    return q

def _validate_terms(terms):
    if not isinstance(terms,dict) or len(terms)>MAX_TERMS: raise ValueError('term budget')
    out={}
    for p,c in terms.items():
        if not isinstance(p,str) or not p or any(x not in 'IXYZ' for x in p): raise ValueError('invalid label')
        q=rat(c)
        if q: out[p]=q
    return out

def _cliques(labels, orders):
    groups=[]; used=set()
    for order in orders:
        for p in order:
            if p in used: continue
            g=[p]
            for q in order:
                if q not in used and q!=p and all(anti(q,z) for z in g): g.append(q)
            groups.append(g); used.update(g)
    # Add overlapping cliques centered on each label.
    for p in labels:
        g=[p]+[q for q in labels if q!=p and anti(p,q)]
        if len(g)>1 and g not in groups: groups.append(g)
    return groups[:MAX_GROUPS]

def propose_fractional_cover(terms, *, groups=None, iterations=64):
    """Return an exact rational overlapping split and bounded construction receipt."""
    t0=perf_counter(); coeff=_validate_terms(terms); labels=list(coeff)
    if groups is None:
        orders=[labels, sorted(labels,key=lambda p:(-abs(float(coeff[p])),p)), list(reversed(labels))]
        groups=_cliques(labels,orders)
    else:
        groups=[list(g) for g in groups]
    if not groups or len(groups)>MAX_GROUPS: raise ValueError('group budget')
    if any(not g for g in groups): raise ValueError('empty group')
    if sum(len(g) for g in groups)>MAX_INCIDENCES: raise ValueError('incidence budget')
    if any(any(p not in coeff for p in g) for g in groups): raise ValueError('unknown term')
    if any(any(not anti(a,b) for a,b in combinations(g,2)) for g in groups): raise ValueError('non-anticommuting group')
    # coordinate descent on floating split, maintaining coefficient sums
    inc={p:[j for j,g in enumerate(groups) if p in g] for p in labels}
    if any(not v for v in inc.values()): raise ValueError('uncovered term')
    x={(j,p):float(coeff[p])/len(inc[p]) for p in labels for j in inc[p]}
    costs=[]
    for it in range(max(0,int(iterations))):
        costs.append(sum(sqrt(sum(x[j,p]**2 for p in g)) for j,g in enumerate(groups)))
        for p in labels:
            js=inc[p]; target=float(coeff[p]); delta=(target-sum(x[j,p] for j in js))/len(js)
            for j in js: x[j,p]+=delta
        # shrink overlap by moving mass from larger norm groups to smaller where possible
        for p in labels:
            js=inc[p]
            if len(js)>1:
                norms={j:sqrt(sum(x[j,q]**2 for q in groups[j])) for j in js}
                jmax=max(js,key=lambda j:norms[j]); jmin=min(js,key=lambda j:norms[j])
                step=0.1*abs(x[jmax,p])
                if x[jmax,p]*x[jmin,p]>=0:
                    x[jmax,p]-=step; x[jmin,p]+=step
    # Exact rationalization: denominator 10^6 then reconcile each term in final incidence.
    split={}
    for j,g in enumerate(groups):
        split[str(j)]={p:Fraction(round(x[j,p]*10**6),10**6) for p in g if round(x[j,p]*10**6)}
    for p in labels:
        js=inc[p]; residual=coeff[p]-sum(split[str(j)].get(p,Fraction(0)) for j in js)
        j=min(js,key=lambda k:len(split[str(k)])); split[str(j)][p]=split[str(j)].get(p,Fraction(0))+residual
        if not split[str(j)][p]: del split[str(j)][p]
    elapsed=perf_counter()-t0
    uppers={str(j):str(Fraction.from_float(sqrt(sum(float(rat(split[str(j)].get(p,0)))**2 for p in g))).limit_denominator(10**9)) for j,g in enumerate(groups)}
    claimed=sum((rat(v) for v in uppers.values()),Fraction(0))
    return {'schema':'v8-fractional-cover-1','terms':{p:str(c) for p,c in coeff.items()},'groups':[g for g in groups],
            'split':split,'upper_bounds':uppers,'claimed_bound':str(claimed),
            'construction_cost':{'iterations':len(costs),'objective_trace':costs,'terms':len(coeff),'incidences':sum(map(len,groups)),'groups':len(groups),'elapsed_seconds':elapsed}}

def _bound(coeff,groups,split):
    total=Fraction(0)
    for j,g in enumerate(groups): total += Fraction.from_float(sqrt(sum(float(rat(split[str(j)].get(p,0)))**2 for p in g))).limit_denominator(10**9)
    return total

def check_fractional_cover(terms,witness):
    t0=perf_counter(); cost={'pair_checks':0,'squares':0,'reconstruction_ops':0}
    try:
        coeff=_validate_terms(terms)
        if not isinstance(witness,dict) or witness.get('schema')!='v8-fractional-cover-1': raise ValueError('schema')
        groups=witness['groups']; split=witness['split']
        if len(coeff)>MAX_TERMS or len(groups)>MAX_GROUPS: raise ValueError('budget')
        if sum(len(g) for g in groups)>MAX_INCIDENCES: raise ValueError('incidence budget')
        if len(split)!=len(groups): raise ValueError('split coverage')
        rec={p:Fraction(0) for p in coeff}
        total=Fraction(0)
        for j,g in enumerate(groups):
            if not isinstance(g,list) or not g: raise ValueError('group')
            ss=split[str(j)]
            for a,b in combinations(g,2):
                cost['pair_checks']+=1
                if not anti(a,b): raise ValueError('non-anticommuting group')
            sq=Fraction(0)
            for p in g:
                if p not in coeff: raise ValueError('unknown label')
                q=rat(ss.get(p,0)); sq+=q*q; rec[p]+=q; cost['reconstruction_ops']+=1
            cost['squares']+=len(g)+1
            # exact rational upper squared bound, allowing decimal/rational upper witness
            upper=rat(witness.get('upper_bounds',{}).get(str(j),str(sqrt(float(sq)))))
            if upper<0 or upper*upper<sq: raise ValueError('invalid upper bound')
            total+=upper
        if rec!=coeff: raise ValueError('exact reconstruction')
        claimed=rat(witness['claimed_bound'])
        if claimed!=total: raise ValueError('claimed mismatch')
        return {'status':'certified','bound':str(total),'checking_cost':cost,'checking_elapsed_seconds':perf_counter()-t0}
    except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError) as e:
        return {'status':'rejected','reason':str(e),'checking_cost':cost,'checking_elapsed_seconds':perf_counter()-t0}
