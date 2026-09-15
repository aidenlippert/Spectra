"""Exact modular injectivity receipt for quartic sector multipliers.

Full column rank modulo a verified prime proves full column rank over Q.
This is an algebraic reduction of equality constraints, not a numerical
singular-value threshold or a lower-bound-family impossibility claim.
"""
from math import isqrt
import time
from experiments.marginal_symbolic import product,number_shift

def rank_mod(rows,columns,prime=65521):
    if prime<2 or any(prime%d==0 for d in range(2,isqrt(prime)+1)):raise ValueError('Prime field required')
    pivots={};used=[]
    for idx,raw in sorted(enumerate(rows),key=lambda p:(len(p[1]),p[0])):
        row={i:int(c)%prime for i,c in raw.items() if int(c)%prime}
        while row:
            i=min(row);a=row[i]
            if i not in pivots:
                inv=pow(a,-1,prime);pivots[i]={j:c*inv%prime for j,c in row.items()};used.append(idx);break
            for j,c in pivots[i].items():
                v=(row.get(j,0)-a*c)%prime
                if v:row[j]=v
                elif j in row:del row[j]
        if len(pivots)==columns:break
    return len(pivots),used

def certify(m,n,basis):
    start=time.monotonic();quartic=[p for p in basis if max(map(len,p),default=0)==4];rows={};prime=65521
    for j,p in enumerate(quartic):
        for w,c in product(number_shift(m,n),p).items():
            if len(w)==6:rows.setdefault(w,{})[j]=(c.numerator*pow(c.denominator,-1,prime))%prime
    ordered=sorted(rows);rank,used=rank_mod([rows[w] for w in ordered],len(quartic),prime)
    if rank!=len(quartic):raise ValueError('Quartic multiplier injectivity not proved')
    return {'kind':'exact_quartic_ideal_injectivity_mod_prime_v1','modes':m,'particles':n,'prime':prime,
            'columns':len(quartic),'rows':len(rows),'rank':rank,'pivot_rows':[list(map(list,ordered[i])) for i in used],
            'seconds':time.monotonic()-start,
            'conclusion':'If a paired degree-three SOS and degree-four H match all coefficients, the quartic part of X must vanish. Drop only these proved zero variables and their now identically zero equations.'}
