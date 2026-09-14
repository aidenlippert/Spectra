"""Independent linear responses to both sides of a cubic anticommutator.

Each abstract PSD block acts through two maps, into charge -1 and +1
polynomials. Cubic coefficients are shared; linear coefficients are free
on each side. Their squares have no sextic contribution. This is a
two-particle-positivity construction related to strengthened T2 conditions.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import time
import numpy as np
from scipy import sparse
from experiments.marginal_symbolic import decode
from experiments.marginal_coefficient import gram_map,dagger
from research.certificate_scaling.adaptive_block_discovery import partition
from research.reconstruction_compression_20260914.inputs import dump,sha
from research.reconstruction_compression_20260914.moments import Proposal

def embed(M,slots,size,transpose=False):
    a=M.tocoo();n=len(slots);i,j=np.divmod(a.col,n)
    if transpose:i,j=j,i
    return sparse.csc_matrix((a.data,(a.row,np.array(slots)[i]*size+np.array(slots)[j])),shape=(M.shape[0],size*size))

def build(case,fixture,state,out):
    t=time.monotonic();folder=Path(out);folder.mkdir(parents=True,exist_ok=False)
    data=json.loads(Path(fixture).read_text());h=decode(data['hamiltonian'],data['modes'],4);m,n=data['modes'],data['particles']
    raw,sig,stats=partition(h,m,'mixed',True)
    raw=[{**g,'words':[w for w in g['words'] if len(w)!=1]} if g['name'].startswith('mixed') else g for g in raw]
    raw=[g for g in raw if g['words']]
    rows=[tuple((1,i) for i in a)+tuple((0,i) for i in b) for k in range(3)
          for a in combinations(range(m),k) for b in combinations(range(m),k)
          if a<=b and sig(tuple((1,i) for i in a)+tuple((0,i) for i in b))==sig(())]
    lookup={w:i for i,w in enumerate(rows)};used=set();groups=[];pairs=[];work=0
    for i,g in enumerate(raw):
        if i in used:continue
        words=g['words'];nc=len(words)
        if max(map(len,words))<=2:
            gid=len(groups);groups.append(g);slots=list(range(nc));pair={'members':[gid],'slot_indices':[slots],'cubic_size':0,'joint_size':nc}
            M=gram_map(words,lookup).tocsc();work+=nc**2
        else:
            found=None
            for j,other in enumerate(raw):
                if j==i or j in used or len(other['words'])!=nc:continue
                loc={w:k for k,w in enumerate(other['words'])}
                if all(dagger(w) in loc for w in words):found=j;break
            if found is None:raise ValueError('Missing adjoint dictionary')
            used.add(found);charge=sum(2*c-1 for c,k in words[0]);signature=sig(words[0])
            linear=[((int(charge>0),k),) for k in range(m) if abs(charge)==1 and sig(((int(charge>0),k),))==signature]
            # Both physical dictionaries use the same cubic index order.
            minus=words+linear;plus=[dagger(w) for w in words]+[dagger(w) for w in linear]
            nl=len(linear);size=nc+2*nl;sminus=list(range(nc+nl));splus=list(range(nc))+list(range(nc+nl,size))
            gid=len(groups);groups.extend([{'name':g['name']+' response','words':minus},{'name':raw[found]['name']+' response','words':plus}])
            pair={'members':[gid,gid+1],'slot_indices':[sminus,splus],'cubic_size':nc,'joint_size':size,'linear_size_per_side':nl}
            M=embed(gram_map(minus,lookup),sminus,size)+embed(gram_map(plus,lookup),splus,size,True);work+=len(minus)**2+len(plus)**2
        M.eliminate_zeros();sparse.save_npz(folder/f'paired_map_{len(pairs)}.npz',M);pairs.append(pair)
    end=time.monotonic();oracle=Proposal(data,json.loads(Path(state).read_text()));oracle.fill(rows)
    weight=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in rows]);y=weight*np.array([oracle.cache[w] for w in rows]);Cs={}
    for k,pair in enumerate(pairs):
        M=sparse.load_npz(folder/f'paired_map_{k}.npz');size=pair['joint_size'];C=np.asarray(M.T@y).reshape(size,size);Cs[f'C_{k}']=(C+C.T)/2
    np.savez_compressed(folder/'moments.npz',**Cs)
    dump(folder/'frame.json',{'kind':'linear_response_v1','case':case,'fixture':str(fixture),'state':str(state),'fixture_sha256':sha(fixture),'state_sha256':sha(state),'groups':groups,'rows':rows,'pairs':pairs,'stats':stats})
    rec={'case':case,'map_seconds':end-t,'moment_seconds':time.monotonic()-end,'total_seconds':time.monotonic()-t,'coefficient_rows':len(rows),
         'requested_moments':len(rows),'maximum_moment_degree':4,'cubic_word_pairs_visited':work,'map_nonzeros':sum(sparse.load_npz(folder/f'paired_map_{k}.npz').nnz for k in range(len(pairs))),
         'minimum_guide_eigenvalue':min(float(np.linalg.eigvalsh(C)[0]) for C in Cs.values()),'many_body_states_enumerated':0,'full_cubic_teacher_used':False}
    dump(folder/'receipt.json',rec);print(json.dumps(rec),flush=True)

def select(meta,Cs,rank):
    maps=[];guide=[]
    for k,p in enumerate(meta['pairs']):
        size=p['joint_size'];nc=p['cubic_size'];C=Cs[f'C_{k}']
        if not nc:J=np.eye(size)
        else:
            nl=size-nc;D=C[nc:,nc:];B=C[:nc,nc:]
            if nl:
                e,U=np.linalg.eigh(D);inv=(U*np.where(e>1e-8,1/np.maximum(e,1e-8),0))@U.T;response=-inv@B.T
            else:response=np.empty((0,nc))
            conditional=C[:nc,:nc]+B@response;conditional=(conditional+conditional.T)/2
            ev,V=np.linalg.eigh(conditional);V=V[:,:min(rank,nc)];r=V.shape[1]
            J=np.zeros((size,r+nl));J[:nc,:r]=V;J[nc:,:r]=response@V;J[nc:,r:]=np.eye(nl)
            J=np.rint(J*1e8)/1e8
            guide.append({'block':k,'cubic_dictionary_size':nc,'retained_cubic_rank':r,'linear_response_size':nl,'smallest_conditional_eigenvalue':float(ev[0]),'last_retained_conditional_eigenvalue':float(ev[r-1]),'largest_response_coefficient':float(abs(response).max()) if nl else 0})
        maps.append(([(p['members'][0],J)],meta['groups'][p['members'][0]]['name']))
    return maps,guide

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('fixture');p.add_argument('state');p.add_argument('out');a=p.parse_args();build(a.case,a.fixture,a.state,a.out)
