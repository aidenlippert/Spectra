"""Exact finite-degree SU(2) projection and an independent row basis."""
from fractions import Fraction as F
from math import lcm,isqrt
import numpy as np
from scipy import sparse
from experiments.marginal_symbolic import mono,adj,canonical,add
from research.certificate_scaling.spin_twirl import twirl

def build(rows):
    lookup={w:i for i,w in enumerate(rows)};ri=[];ci=[];values=[]
    for j,w in enumerate(rows):
        a=canonical(adj(mono(w)));p=mono(w) if a==mono(w) else add(mono(w),a)
        q=twirl(p)
        for v,c in q.items():
            if v in lookup:ri.append(lookup[v]);ci.append(j);values.append(F(c))
            elif not all(z in lookup for z in canonical(adj(mono(v)))):raise ValueError('Spin projector leaves declared row space')
    den=lcm(*(c.denominator for c in values));iv=np.array([int(c*den) for c in values],dtype=np.int64)
    A=sparse.csr_matrix((iv,(ri,ci)),shape=(len(rows),len(rows)))
    bound=int(abs(iv).max())**2*len(rows)
    if bound>=2**63:raise ValueError('Integer projector check would overflow')
    defect=A@A-den*A;defect.eliminate_zeros()
    if defect.nnz:raise AssertionError('Exact spin projector is not idempotent')
    trace=F(int(A.diagonal().sum()),den)
    if trace.denominator!=1:raise AssertionError('Projector trace is not integral')
    # Independence modulo a prime implies independence over Q. Since an
    # idempotent has rank equal to its trace, trace-many independent rows
    # constitute a rational row basis; none of the discarded equations is lost.
    prime=65521
    if any(prime%d==0 for d in range(2,isqrt(prime)+1)):raise AssertionError('Rank modulus is not prime')
    pivots={};selected=[]
    for i in range(A.shape[0]):
        row=A.getrow(i);d={int(j):int(v)%prime for j,v in zip(row.indices,row.data) if int(v)%prime}
        while d:
            p=min(d)
            if p not in pivots:
                inv=pow(d[p],prime-2,prime);pivots[p]={j:v*inv%prime for j,v in d.items()};selected.append(i);break
            factor=d[p]
            for j,v in pivots[p].items():
                x=(d.get(j,0)-factor*v)%prime
                if x:d[j]=x
                else:d.pop(j,None)
    if len(selected)!=trace:raise AssertionError('Modular row rank differs from exact projector trace')
    T=A.astype(float)/den
    return T,np.array(selected),{'original_rows':len(rows),'independent_invariant_rows':len(selected),'integer_denominator':den,
        'projector_nonzeros':A.nnz,'exact_integer_idempotence_checked':True,'modular_prime':prime,'rank_equals_exact_trace':True}
