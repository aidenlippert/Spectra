"""Exact maximum charge-metric ratios by finite-memory dynamic programming."""
from fractions import Fraction as F


def ratio_upper(oracle, delta, mask, bits, minimum_doublons=1):
    """Maximize one ratio over a source event and a fixed-spin ionic tail.

    Charge factors must have one- or two-site support. Positive products are
    maximized exactly; no occupation configurations are listed or acted on.
    None denotes an empty source domain. Costs count all visited DP states.
    """
    m=oracle.sites
    if (type(minimum_doublons) is not int or not 0<=minimum_doublons<=m//2 or
        type(mask) is not int or type(bits) is not int or not 0<=bits<=mask<=oracle.all_bits or bits&~mask or
        type(delta) not in (tuple,list) or len(delta)!=m or any(type(d) is not int for d in delta)):
        raise ValueError('Valid charge change, occupation event and doublon threshold required')
    factors=[(label,value) for label,value in oracle.local_factors if value!=1 and any(delta[i] for i,_ in label)]
    if any(not 1<=len(label)<=2 for label,_ in factors):
        raise ValueError('One- and two-site charge factors required')
    memory=max((max(i for i,_ in label)-min(i for i,_ in label) for label,_ in factors),default=0)
    if memory>4:raise ValueError('Exact ratio DP currently limited to charge range four')
    ends=[[] for _ in range(m)]
    for label,value in factors:ends[max(i for i,_ in label)].append((label,value))
    dp={(0,0,0,()):F(1)};visited=1;transitions=0;peak=1;scalar_powers=0;cache={}
    for site in range(m):
        local_mask=(mask>>(2*site))&3;local_bits=(bits>>(2*site))&3
        options=[value for value in range(4) if value&local_mask==local_bits]
        nxt={}
        for (a,b,d,tail),weight in dp.items():
            for value in options:
                aa=a+(value&1);bb=b+((value>>1)&1)
                remaining=m-site-1
                if aa>oracle.target or bb>oracle.target or aa+remaining<oracle.target or bb+remaining<oracle.target:continue
                dd=min(minimum_doublons,d+(value==3))
                if dd+remaining<minimum_doublons:continue
                charge=value.bit_count()-1;history=tail+(charge,)
                key=(site,history)
                if key not in cache:
                    factor=F(1)
                    for label,base in ends[site]:
                        before=after=1
                        for i,power in label:
                            q=history[-1-(site-i)]
                            before*=q**power;after*=(q+delta[i])**power
                        factor*=base**(after-before);scalar_powers+=1
                    cache[key]=factor
                candidate=weight*cache[key]
                state=(aa,bb,dd,history[-memory:] if memory else ())
                if state not in nxt or candidate>nxt[state]:nxt[state]=candidate
                transitions+=1
        dp=nxt;visited+=len(dp);peak=max(peak,len(dp))
    values=[v for (a,b,d,_),v in dp.items() if a==b==oracle.target and d==minimum_doublons]
    return (max(values) if values else None), {'visited_states':visited,'transitions':transitions,
        'peak_states':peak,'memory_sites':memory,'metric_scalar_powers':scalar_powers,
        'occupation_endpoint_evaluations':0}
