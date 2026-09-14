"""Exact charge-pattern transition matrix and PF obstruction diagnostic."""
from itertools import product
from fractions import Fraction
import json
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigs

def patterns(sites=8):
    return [q for q in product((-1,0,1), repeat=sites) if sum(q)==0 and any(q)]

def build(sites=8):
    ps=patterns(sites); index={q:i for i,q in enumerate(ps)}
    A=lil_matrix((len(ps),len(ps)),dtype=float)
    for r,q in enumerate(ps):
        d=sum(x*x for x in q)//2
        # Generator-style comparison operator: diagonal penalty, positive
        # transition adjacency.  Its Perron root is the sharp charge-only
        # threshold for inequalities A w >= 0.
        A[r,r]=-4*d
        for i in range(sites-1):
            a,b=q[i]+1,q[i+1]+1
            targets=[]
            if abs(a-b)==1:
                z=list(q);z[i],z[i+1]=z[i+1],z[i];targets.append((tuple(z),1))
            elif {a,b}=={0,2}:
                z=list(q);z[i]=z[i+1]=0;targets.append((tuple(z),2))
            elif a==b==1:
                for direction in (-1,1):
                    z=list(q);z[i]=direction;z[i+1]=-direction;targets.append((tuple(z),1))
            for z,w in targets:
                if z in index:A[r,index[z]] += w
    return ps,A.tocsr()

def run(gammas=(-18,-16,-14,-6),sites=8):
    ps,A=build(sites); vals,vecs=eigs(A.T,k=1,which='LR'); lam=float(vals[0].real)
    y=np.abs(vecs[:,0].real); y/=y.sum()
    scale=10**12; yi=np.maximum(1,np.rint(y*scale).astype(object));
    out=[]
    for gamma in gammas:
        z=np.asarray(yi,dtype=object) @ (A.toarray().astype(object)-gamma*np.eye(len(ps),dtype=object))
        out.append({'gamma':gamma,'max_component':int(max(z)),'min_component':int(min(z)),'strict_negative':bool(max(z)<0)})
    result={'sites':sites,'patterns':len(ps),'nnz':int(A.nnz),'pf_eigenvalue':lam,'rounded_scale':scale,'tests':out}
    p=Path('results/marginal_graded_hubbard8/charge_only_limit');p.mkdir(parents=True,exist_ok=True)
    (p/'diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__': print(json.dumps(run(),indent=2))
