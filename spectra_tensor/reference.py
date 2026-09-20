"""Matched-accuracy full-sector control, with memory preflight before allocation.

Uses a Kronecker-sum hopping action, avoiding a full many-body Hamiltonian.
Its solution vectors still enumerate the complete balanced-spin sector.
"""
from itertools import combinations
from math import comb
from time import perf_counter
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import gmres,LinearOperator
from .exact import validate_model

def solve_reference(spec,z,target=.001,max_bytes=1_000_000_000,restart=32,maxiter=20):
    start=perf_counter();n,U,edges=validate_model(spec);m=comb(n,n//2);dimension=m*m
    if target<=0 or z.imag<=0:raise ValueError('positive target/broadening')
    estimate=(restart+12)*16*dimension+8*dimension
    if estimate>max_bytes:return dict(status='resource_limit_preflight',dimension=dimension,single_complex_vector_bytes=16*dimension,projected_working_bytes=estimate,max_bytes=max_bytes,seconds=perf_counter()-start)
    labels=sorted(sum(1<<i for i in c) for c in combinations(range(n),n//2));lookup={s:i for i,s in enumerate(labels)};rows=[];cols=[];values=[]
    for j,state in enumerate(labels):
        for a,b,t in edges:
            if ((state>>a)&1)==((state>>b)&1):continue
            parity=(state&((1<<b)-(1<<(a+1)))).bit_count()%2
            rows.append(lookup[state^(1<<a)^(1<<b)]);cols.append(j);values.append(-float(t)*((-1)**parity))
    hopping=csr_matrix((values,(rows,cols)),shape=(m,m));diag=np.zeros((m,m))
    for i,u in enumerate(U):
        occ=np.array([(s>>i)&1 for s in labels],float);diag+=float(u)*occ[:,None]*occ[None,:]
    calls=[0]
    def action(v):
        X=v.reshape(m,m);calls[0]+=1;return (hopping@X+(hopping@X.T).T+diag*X).reshape(-1)
    b=np.zeros(dimension,complex);up=sum(1<<j for j in range(0,n,2));down=sum(1<<j for j in range(1,n,2));index=lookup[up]*m+lookup[down];b[index]=1
    A=LinearOperator((dimension,dimension),matvec=lambda v:z*v-action(v),dtype=complex);pre=z-diag.reshape(-1);M=LinearOperator((dimension,dimension),matvec=lambda v:v/pre,dtype=complex)
    x,info=gmres(A,b,M=M,rtol=.5*np.sqrt(target*z.imag),atol=0.,restart=restart,maxiter=maxiter)
    hx=action(x);r=b-z*x+hx;center=2*x[index]-x@(z*x-hx);radius=float(np.vdot(r,r).real/z.imag)
    return dict(status='target_met_numerically' if radius<=target else 'target_not_met',info=int(info),dimension=dimension,center_real=float(center.real),center_imag=float(center.imag),radius=radius,matvecs=calls[0],seconds=perf_counter()-start,projected_working_bytes=estimate,acceptance='numerical residual only; no independent exact receipt')
