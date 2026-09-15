"""Share spin-flipped factors inside an exactly spin-averaged proof."""
import numpy as np
from experiments.marginal_symbolic import canonical,mono
from experiments.marginal_coefficient import dagger

def locations(words):
    out={}
    for i,w in enumerate(words):
        p=canonical(mono(w))
        if len(p)!=1:raise ValueError('Monomial dictionary required')
        cw,s=next(iter(p.items()));out[cw]=(i,int(s))
    return out

def flipped(word):
    p=canonical(mono(tuple((c,i^1) for c,i in word)))
    if len(p)!=1:raise ValueError('Spin flip did not remain a monomial')
    return next(iter(p.items()))

def member_maps(groups,old,V):
    i=old[0][0];members=[(i,V)]
    if len(old)==2:
        j=old[1][0];loc={w:k for k,w in enumerate(groups[j]['words'])};order=[loc[dagger(w)] for w in groups[i]['words']]
        W=np.empty_like(V);W[order]=V;members.append((j,W))
    return members

def share(groups,maps,tolerance=1e-6):
    locs=[locations(groups[mem[0][0]]['words']) for mem,name in maps]
    output=[];ids=[];used=set();diagnostics=[]
    for k,(members,name) in enumerate(maps):
        if k in used:continue
        words=groups[members[0][0]]['words'];V=members[0][1];target={flipped(w)[0] for w in words}
        matches=[j for j,loc in enumerate(locs) if set(loc)==target and len(maps[j][0])==len(members)]
        if len(matches)!=1:raise ValueError(('Spin-flip orbit is ambiguous',name,matches))
        partner=matches[0]
        if partner==k:
            n=V.shape[0]
            if not np.array_equal(V,np.eye(n)):raise ValueError('Self-flip splitting requires the full coordinate block')
            used_rows=set();sectors={-1:[],1:[]}
            for a,w in enumerate(words):
                if a in used_rows:continue
                cw,s=flipped(w);b,t=locs[k][cw];sign=int(s/t)
                if a==b:
                    c=np.zeros(n);c[a]=1;sectors[sign].append(c)
                else:
                    cw2,s2=flipped(words[b]);back,t2=locs[k][cw2]
                    if back!=a or int(s2/t2)!=sign:raise AssertionError('Spin flip is not an involution')
                    for parity in (-1,1):
                        c=np.zeros(n);c[a]=0.70710678;c[b]=parity*sign*0.70710678;sectors[parity].append(c)
                used_rows.update((a,b))
            for parity,cols in sectors.items():
                if cols:
                    B=np.column_stack(cols);output.append((member_maps(groups,members,B),name+f' spin parity {parity}'));ids.append(k)
            diagnostics.append({'block':k,'self_flip':True,'before_dimension':n,'after_dimensions':[len(v) for v in sectors.values() if v]})
        else:
            other=maps[partner][0][0];other_words=groups[other[0]]['words'];W=np.empty((len(words),other[1].shape[1]))
            for a,w in enumerate(other_words):
                cw,s=flipped(w);b,t=locs[k][cw];W[b]=int(s/t)*other[1][a]
            if np.array_equal(V,np.eye(V.shape[0])) and V.shape[1]==W.shape[1]:B=V;singular=[]
            else:
                U,s,_=np.linalg.svd(np.column_stack((V,W)),full_matrices=False);keep=s>tolerance;B=np.rint(U[:,keep]*1e8)/1e8;singular=s.tolist()
            output.append((member_maps(groups,members,B),name+' shared spin orbit'));ids.append(k);used.add(partner)
            diagnostics.append({'block':k,'partner':partner,'before_dimensions':[V.shape[1],W.shape[1]],'after_dimension':B.shape[1],
                'span_tolerance_proposal_only':tolerance,'union_singular_values':singular})
    return output,ids,diagnostics
