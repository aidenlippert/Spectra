"""Column generation of coherent positive squares in a quartic paired cone.

Quadratic PSD blocks stay full. Cubic blocks use nonnegative combinations of
fixed polynomial factors, and full coefficient-dual violations add new factors.
No many-body state or matrix is formed.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from experiments.marginal_coefficient import dagger
from research.collective_completion_20260914.linear_response import embed_paired_seed


def paired(groups,members,V):
    i=members[0][0];j=members[1][0]
    loc={w:k for k,w in enumerate(groups[j]['words'])};order=[loc[dagger(w)] for w in groups[i]['words']]
    W=np.empty_like(V);W[order]=V
    return [(i,V),(j,W)]


def from_guide(groups,maps):
    result=[];ids=[];receipt=[]
    for k,(members,name) in enumerate(maps):
        if len(members)==1:result.append((members,name));ids.append(k);continue
        V=members[0][1]
        for j in range(V.shape[1]):
            result.append((paired(groups,members,V[:,j:j+1]),name+f' guide atom {j}'));ids.append(k)
        receipt.append({'source_map':k,'guide_atoms':V.shape[1]})
    return result,ids,receipt


def initial(meta,groups,maps,seed):
    path=Path(seed);d=json.loads((path.parent/'construction.json').read_text())
    old_meta=json.loads((Path(d['prepared_dependency'])/'frame.json').read_text())
    if old_meta.get('kind')!='paired_quartic_v1' or d.get('shared_spin_orbits'):raise ValueError('Unshared paired seed required')
    if old_meta['rows']!=meta['rows'] or old_meta['fixture_sha256']!=meta['fixture_sha256']:raise ValueError('Different coefficient problem')
    cert=json.loads((path.parent/'certificate.json').read_text())['core'];den=cert['denominator']
    result=[];ids=[];receipt=[]
    for k,(members,name) in enumerate(maps):
        if len(members)==1:result.append((members,name));ids.append(k);continue
        old=old_meta['groups'][old_meta['pairs'][k]['members'][0]]
        words=[tuple(tuple(x) for x in w) for w in old['words']]
        signature=frozenset(words)
        matches=[b for b in cert['blocks'] if frozenset(tuple(tuple(x) for x in w) for w in b['words'])==signature]
        if len(matches)>1:raise ValueError(('Ambiguous certified seed block',old['name']))
        if not matches:
            v=np.zeros((len(groups[members[0][0]]['words']),1));v[0,0]=1
            result.append((paired(groups,members,v),name+' initially unused atom'));ids.append(k)
            receipt.append({'source_map':k,'source_rows':0,'old_point_weights':[0.],
                            'reason':'This dictionary contributes no factor in the certified seed'})
            continue
        b=matches[0];loc={tuple(tuple(x) for x in w):i for i,w in enumerate(b['words'])}
        L=np.array(b['factor'],dtype=float)[:,[loc[w] for w in words]]/den
        new_words=groups[members[0][0]]['words'];size=len(new_words)
        V=embed_paired_seed(words,L.T,new_words,{'joint_size':size,'cubic_size':0,'slot_indices':[list(range(size))]})
        weights=[]
        for col in V.T:
            norm=float(np.linalg.norm(col))
            if not norm:continue
            # Power-of-two rescaling changes conditioning, not the intended ray.
            unit=2.**np.floor(np.log2(norm));v=(col/unit)[:,None]
            result.append((paired(groups,members,v),name+f' certified atom {len(weights)}'));ids.append(k);weights.append(float(unit**2))
        receipt.append({'source_map':k,'source_rows':L.shape[0],'old_point_weights':weights,
                        'old_dictionary_size':len(words),'new_dictionary_size':size})
    return result,ids,receipt


def add(meta,groups,maps,ids,prepared,y,count,max_atoms=4000):
    result=list(maps);out_ids=list(ids);details=[]
    for k,p in enumerate(meta['pairs']):
        if len(p['members'])!=2:continue
        first=next(mem for (mem,name),mid in zip(maps,ids) if mid==k)
        M=sparse.load_npz(Path(prepared)/f'paired_map_{k}.npz');n=first[0][1].shape[0]
        C=np.asarray(M.T@y).reshape(n,n);C=(C+C.T)/2;ev,V=np.linalg.eigh(C)
        take=np.where(ev < -1e-7)[0][:count]
        for j in take:
            v=np.rint(V[:,j:j+1]*1e8)/1e8
            result.append((paired(groups,first,v),f'map {k} dual atom {len(result)}'));out_ids.append(k)
        details.append({'map':k,'lowest_full_dual_eigenvalue':float(ev[0]),'added_atoms':len(take)})
    if len(result)>max_atoms:raise ValueError('Atom-count budget exceeded')
    return result,out_ids,details
