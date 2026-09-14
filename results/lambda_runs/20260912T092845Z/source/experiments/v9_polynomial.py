"""Conventional projected ODE collocation and matched projected Taylor.

Floating linear algebra proposes rational polynomials. No spectral error
estimate is trusted; the original Pauli residual checker supplies the bound.
"""
from fractions import Fraction as F
from time import perf_counter
import numpy as np
from .v7_certificate import (Piece,clean,rational,BudgetExceeded,residual_records,norm_witness)
GRID=10**12


def quantize(value):
    if not np.isfinite(value):raise ValueError('nonfinite polynomial proposal')
    return rational(F(int(round(float(value)*GRID)),GRID))


def krylov_basis(gen,initial,dimension=8):
    if type(dimension) is not int or not 1<=dimension<=32:
        raise ValueError('Krylov dimension cap')
    initial=clean(initial,gen.n,gen.max_terms)
    beta=sum(float(c)**2 for c in initial.values())**.5
    if not beta or not np.isfinite(beta):raise ValueError('nonzero finite initial required')
    vectors=[{p:float(c)/beta for p,c in initial.items()}]
    columns=[];orth=0;conversions=0
    for j in range(dimension):
        q={p:quantize(c) for p,c in vectors[j].items()};conversions+=len(q)
        w={p:float(c) for p,c in gen.apply(q).items()}
        h=[0.]*(j+2)
        for _ in range(2):
            for i,v in enumerate(vectors):
                dot=sum(w.get(p,0.)*c for p,c in v.items());orth+=len(v)
                h[i]+=dot
                for p,c in v.items():w[p]=w.get(p,0.)-dot*c;orth+=1
        norm=sum(c*c for c in w.values())**.5;orth+=len(w)
        if not np.isfinite(norm):raise ValueError('nonfinite Krylov proposal')
        h[j+1]=norm;columns.append(h)
        if norm<=1e-13 or j+1==dimension:break
        vectors.append({p:c/norm for p,c in w.items() if c})
    k=len(vectors);labels=sorted(set().union(*(v.keys() for v in vectors)))
    if len(labels)>gen.max_terms:raise BudgetExceeded('Krylov support cap')
    V=np.array([[v.get(p,0.) for v in vectors] for p in labels],dtype=float)
    A=np.array([[columns[j][i] if i<len(columns[j]) else 0. for j in range(k)] for i in range(k)],dtype=float)
    return dict(labels=labels,V=V,A=A,beta=beta,cost=dict(orthogonalization_scalar_ops=orth,quantized_entries=conversions,dimension=k,stored_float_entries=V.size+A.size,generator=dict(gen.cost)))


def _piece(basis,initial,duration,normalized_coefficients):
    T=rational(duration)
    if T<=0:raise ValueError('duration')
    cs=[dict(initial)]
    for ell,c in enumerate(normalized_coefficients,1):
        op={p:quantize(x)/T**ell for p,x in zip(basis['labels'],basis['V']@c)}
        cs.append({p:v for p,v in op.items() if v})
    if sum(map(len,cs))>50000:raise BudgetExceeded('polynomial coefficient budget')
    return Piece(T,tuple(cs))


def projected_taylor_piece(basis,initial,duration,degree):
    if type(degree) is not int or not 0<=degree<=24:raise ValueError('degree cap')
    T=rational(duration);A=basis['A'];k=A.shape[0]
    c=np.zeros(k);c[0]=basis['beta'];coeff=[]
    for ell in range(1,degree+1):
        c=(float(T)/ell)*(A@c);coeff.append(c.copy())
    return _piece(basis,initial,T,coeff),dict(projected_multiply_adds=degree*k*k,reconstruction_multiply_adds=degree*basis['V'].size,linear_system_dimension=0)


def collocation_piece(basis,initial,duration,degree):
    """Gauss-node derivative collocation in a fixed Krylov subspace.

U_i=y0+T sum_j Q_ij A U_j; P(T*s)=y0+T sum_j int_0^s l_j(v)dv A U_j.
All node/solve errors are included in the exported polynomial's residual.
"""
    if type(degree) is not int or not 1<=degree<=12:raise ValueError('collocation degree cap')
    T=rational(duration)
    if T<=0:raise ValueError('duration')
    A=basis['A'];k=A.shape[0];m=degree
    if k*m>384:raise BudgetExceeded('collocation dense solve dimension')
    nodes=(np.polynomial.legendre.leggauss(m)[0]+1.)/2.
    L=np.zeros((m,m));Q=np.zeros((m,m))
    for j,node in enumerate(nodes):
        coeff=np.array([1.]);denom=1.
        for ell,other in enumerate(nodes):
            if ell!=j:
                coeff=np.polynomial.polynomial.polymul(coeff,np.array([-other,1.]));denom*=node-other
        coeff=coeff/denom;L[:,j]=coeff
        integral=np.r_[0.,coeff/np.arange(1,m+1)]
        Q[:,j]=np.polynomial.polynomial.polyval(nodes,integral)
    system=np.eye(m*k)-float(T)*np.kron(Q,A)
    y0=np.zeros(k);y0[0]=basis['beta']
    stages=np.linalg.solve(system,np.tile(y0,m)).reshape(m,k)
    velocities=A@stages.T
    cs=[(float(T)/ell)*(velocities@L[ell-1,:]) for ell in range(1,m+1)]
    return _piece(basis,initial,T,cs),dict(linear_system_dimension=m*k,linear_system_entries=system.size,
                dense_solve_cubic_size=(m*k)**3,projected_multiply_adds=m*k*k+m*m*k,
                reconstruction_multiply_adds=m*basis['V'].size,stage_float_entries=stages.size)


def best_witness(gen,initial,piece,tolerance,integration_basis='bernstein'):
    """Reuse recomputed residual records across ordinary norm proposals."""
    start=perf_counter();tol=rational(tolerance)
    records=residual_records(gen,initial,[piece],integration_basis)
    best=None;probes=[];work=dict(group_comparisons=0,sorted_terms=0,sqrt_enclosures=0)
    for mode in ('l1','firstfit','weighted'):
        bound=F(0);ws={}
        for label,op,weight in records:
            groups,cost=norm_witness(op,mode);ws[label]=groups
            bound+=weight*sum((F(g['upper']) for g in groups),F(0))
            for key in work:work[key]+=cost[key]
        probes.append(dict(grouping=mode,bound=str(bound)))
        w=dict(schema='v7-residual-1',integration_basis=integration_basis,witnesses=ws,claimed_bound=str(bound))
        if best is None or bound<F(best['claimed_bound']):best=w
        if bound<=tol:break
    best['construction_cost']=work
    return best,dict(probes=probes,generator=dict(gen.cost),norm_work=work,witness_seconds=perf_counter()-start)
