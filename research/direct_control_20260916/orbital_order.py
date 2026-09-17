"""Direct fermionic orbital permutation of the compact input MPS."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from .validated_tensor import *


def group_spatial(ts):
    out=[];error=0.0
    pref,suf=tensor_bounds(ts)
    for i in range(0,len(ts),2):
        a,b=ts[i:i+2]
        z,e=mm(a.reshape(-1,a.shape[2]),b.reshape(b.shape[0],-1))
        out.append(z.reshape(a.shape[0],4,b.shape[2]))
        error=plus(error,times(pref[i],e,suf[i+2]))
    return out,error


def swap_neighbors(ts,i):
    pref=prefix_gram_bounds(ts);_,suf=tensor_bounds(ts)
    a,b=ts[i:i+2];l,_,_=a.shape;r=b.shape[2]
    z,e=mm(a.reshape(l*4,-1),b.reshape(b.shape[0],4*r))
    z=z.reshape(l,4,4,r)
    parity=[0,1,1,0]
    z=z.transpose(0,2,1,3).copy()
    for s in range(4):
        for t in range(4):
            if parity[s]*parity[t]:z[:,s,t,:]*=-1
    z=z.reshape(l*4,4*r)
    u,ss,vh=svd(z,full_matrices=False,check_finite=False)
    # Keep every direction. Compression is separately charged after ordering.
    rmat=u.conj().T@z
    rec,re=mm(u,rmat)
    de=plus(e,difference_norm(z,rec),re)
    result=list(ts);result[i]=u.reshape(l,4,-1);result[i+1]=rmat.reshape(u.shape[1],4,r)
    return result,times(pref[i],de,suf[i+2])


def permute(ts,order):
    current=list(range(len(ts)));out=ts;error=0.;swaps=0
    for target,orb in enumerate(order):
        at=current.index(orb)
        while at>target:
            out,e=swap_neighbors(out,at-1);error=plus(error,e)
            current[at-1],current[at]=current[at],current[at-1];at-=1;swaps+=1
    return out,error,swaps


def split_spatial(ts):
    pref=prefix_gram_bounds(ts);_,suf=tensor_bounds(ts);out=[];error=0.
    for i,a in enumerate(ts):
        l,p,r=a.shape;z=a.reshape(l*2,2*r)
        u,s,vh=svd(z,full_matrices=False,check_finite=False)
        v=u.conj().T@z;rec,re=mm(u,v)
        e=plus(difference_norm(z,rec),re)
        error=plus(error,times(pref[i],e,suf[i+1]))
        out.extend([u.reshape(l,2,-1),v.reshape(u.shape[1],2,r)])
    return out,error


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--order',type=int,nargs='+',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    outdir=Path(a.out);outdir.mkdir(parents=True,exist_ok=False)
    base=Path('results/direct_control_20260916/imported/Spectra_control_reduction/inputs')
    start=time.monotonic();cert=json.loads((base/'state.json').read_text());ts,e0=load_mps(cert)
    ts,e1=right_canonicalize(ts);ts,e2=group_spatial(ts);ts,e3,swaps=permute(ts,a.order)
    ts,e4=split_spatial(ts);ts,e5=right_canonicalize(ts)
    error=plus(e0,e1,e2,e3,e4,e5)
    np.savez_compressed(outdir/'state.npz',**{f'a{i}':x for i,x in enumerate(ts)})
    record={'order':a.order,'swaps':swaps,'source_permutation_error_bound':error,'preparation_seconds':time.monotonic()-start,'compressions':[]}
    for bond in [16,24,32,48,64,96,128,192]:
        t=time.monotonic();bs,e,loc=compress(ts,bond)
        rr={'bond':bond,'additional_error_bound':e,'total_error_bound':plus(e,error),'seconds':time.monotonic()-t,'norm_interval':expectation(bs)}
        record['compressions'].append(rr);print(json.dumps(rr),flush=True)
    record['total_seconds']=time.monotonic()-start
    (outdir/'record.json').write_text(json.dumps(record,indent=2)+'\n')
