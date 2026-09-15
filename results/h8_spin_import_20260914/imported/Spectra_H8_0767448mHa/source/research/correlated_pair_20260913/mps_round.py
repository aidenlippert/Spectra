"""Numerical charge gauge proposal; exact acceptance lives in mps_exact.py."""
import argparse
import json
from pathlib import Path
import numpy as np
from research.molecular_collective_20260913.core import digest


def rationalize(data,tensors,spin_counts,denominator=10**9):
    A=[np.asarray(x,dtype=float).copy() for x in tensors];m=data['modes']
    if len(A)!=m or any(x.ndim!=3 or x.shape[1]!=2 for x in A): raise ValueError('left,physical,right tensor layout')
    if A[0].shape[0]!=1 or A[-1].shape[2]!=1: raise ValueError('Open MPS required')
    # Left canonicalization uses only local tensor/bond matrices.
    for i in range(m-1):
        l,_,r=A[i].shape; Q,R=np.linalg.qr(A[i].reshape(l*2,r),mode='reduced')
        A[i]=Q.reshape(l,2,Q.shape[1]);A[i+1]=np.einsum('ab,bsc->asc',R,A[i+1])
    nrm=float(np.linalg.norm(A[-1]));A[-1]/=nrm
    base=m+1;charges=[[[0,0]]];worst=0.;discarded=0.
    for i in range(m):
        q=np.array([base*a+b for a,b in charges[-1]],float)
        inc=base if i%2==0 else 1
        if i<m-1:
            T=sum(A[i][:,s,:].T@((q+s*inc)[:,None]*A[i][:,s,:]) for s in (0,1))
            vals,R=np.linalg.eigh((T+T.T)/2);A[i]=np.einsum('asb,bc->asc',A[i],R)
            A[i+1]=np.einsum('ab,bsc->asc',R.T,A[i+1])
            labels=np.rint(vals).astype(int);worst=max(worst,float(np.max(np.abs(vals-labels))))
            right=[[int(v//base),int(v%base)] for v in labels]
        else:right=[list(spin_counts)]
        for a,left in enumerate(charges[-1]):
            for s in (0,1):
                target=list(left);target[i%2]+=s
                for b,rq in enumerate(right):
                    if target!=rq or any(v<0 or v>spin_counts[k] for k,v in enumerate(rq)):
                        discarded+=float(A[i][a,s,b]**2);A[i][a,s,b]=0.
        charges.append(right)
    edges=[]
    for tensor in A:
        rounded=np.rint(tensor*denominator).astype(np.int64)
        edges.append([[int(a),int(s),int(b),int(rounded[a,s,b])] for a,s,b in zip(*np.nonzero(rounded))])
    # Remove all bond vertices outside a nonzero complete boundary path.
    forward=[{0}]
    for e in edges:forward.append({b for a,s,b,v in e if a in forward[-1]})
    backward=[set() for _ in range(m+1)];backward[-1]={0}
    for i in reversed(range(m)):backward[i]={a for a,s,b,v in edges[i] if b in backward[i+1]}
    keep=[sorted(f&b) for f,b in zip(forward,backward)]
    if any(not v for v in keep): raise ValueError('Charge truncation produced zero state')
    maps=[{v:j for j,v in enumerate(k)} for k in keep]
    cert={'kind':'integer_charge_mps_v1','fixture_sha256':digest(data),'modes':m,'particles':data['particles'],
          'spin_counts':list(spin_counts),'denominator':denominator,
          'bond_charges':[[charges[i][a] for a in k] for i,k in enumerate(keep)],
          'tensors':[[[maps[i][a],s,maps[i+1][b],v] for a,s,b,v in e if a in maps[i] and b in maps[i+1]] for i,e in enumerate(edges)]}
    return cert,{'charge_eigenvalue_max_distance_to_integer':worst,'discarded_local_tensor_squared_entries':discarded,
                 'input_norm_before_normalization':nrm,'postround_bond_dimensions':list(map(len,keep)),
                 'warning':'Discarded local tensor entries are a proposal diagnostic, not a state norm or energy error guarantee.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture');p.add_argument('proposal');p.add_argument('output');p.add_argument('--alpha',type=int,required=True);p.add_argument('--beta',type=int,required=True);a=p.parse_args()
    data=json.loads(Path(a.fixture).read_text());proposal=json.loads(Path(a.proposal).read_text())
    cert,stats=rationalize(data,proposal['tensors'],[a.alpha,a.beta]);Path(a.output).write_text(json.dumps(cert,separators=(',',':'))+'\n')
    Path(a.output+'.rounding.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats))
