"""Sparse coefficient SDP assembly exposing both certificate and moment proposals."""
from fractions import Fraction as F
from itertools import combinations
import time

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

from experiments.marginal_coefficient import gram_map
from experiments.marginal_symbolic import canonical, hermitian, multiplier_basis, number_shift, product, validate_word


def assemble_and_solve(h,modes,particles,blocks,degree=3,groups=None,eps=1e-9,residual_penalty=False,parity_masks=None):
    if degree not in (3,4) or type(degree) is not int:
        raise ValueError('Only degree three and four are supported')
    if type(modes) is not int or modes<1 or type(particles) is not int or not 0<=particles<=modes:
        raise ValueError('Invalid sector')
    h=canonical(h)
    if not hermitian(h) or any(len(w)>4 or sum(2*c-1 for c,_ in w) for w in h):
        raise ValueError('Expected Hermitian two-body Hamiltonian')
    for w in h:validate_word(w,modes,4)
    if groups is None:groups=[list(range(modes))]
    if sorted(i for group in groups for i in group)!=list(range(modes)):
        raise ValueError('Groups must partition the modes')
    parity_masks=[] if parity_masks is None else list(parity_masks)
    if any(type(s) is not int or not 0<s<(1<<modes) for s in parity_masks):
        raise ValueError('Invalid parity symmetry mask')
    member={i:g for g,group in enumerate(groups) for i in group}
    def charge(w):
        result=[0]*len(groups)
        for c,i in w:result[member[i]]+=2*c-1
        return tuple(result)
    parity=lambda w:tuple(sum((s>>i)&1 for c,i in w)%2 for s in parity_masks)
    invariant=lambda w:not any(charge(w)) and not any(parity(w))
    if any(not invariant(w) for w in h):raise ValueError('Hamiltonian breaks supplied charges')
    start=time.monotonic();rows=[]
    for k in range(degree+1):
        sets=list(combinations(range(modes),k))
        for i,left in enumerate(sets):
            for right in sets[i:]:
                w=tuple((1,p) for p in left)+tuple((0,p) for p in right)
                if invariant(w):rows.append(w)
    lookup={w:i for i,w in enumerate(rows)}
    basis=[p for p in multiplier_basis(modes,degree-1) if all(invariant(w) for w in p)]
    ri=[];ci=[];values=[]
    shift=number_shift(modes,particles)
    for j,p in enumerate(basis):
        for w,c in product(shift,p).items():
            if w in lookup:ri.append(lookup[w]);ci.append(j);values.append(float(c))
    a=sp.csr_matrix((values,(ri,ci)),shape=(len(rows),len(basis)))
    b=cp.Variable();x=cp.Variable(len(basis));unit=np.zeros(len(rows));unit[lookup[()]]=1
    expression=b*unit+a@x;variables=[];maps=[];dimensions=[];transforms=[];solver_maps=[]
    for block in blocks:
        words=[validate_word(w,modes,degree) for w in block['words']]
        if not words or len({charge(w)+parity(w) for w in words})!=1:
            raise ValueError('Each Gram block must have a single conserved charge vector')
        transform=block.get('basis_transform')
        if transform is not None:
            transform=np.asarray(transform,dtype=float)
            if transform.ndim!=2 or transform.shape[0]!=len(words) or not transform.shape[1] or not np.all(np.isfinite(transform)):
                raise ValueError('Invalid polynomial basis transform')
        dimension=len(words) if transform is None else transform.shape[1]
        g=gram_map(words,lookup);q=cp.Variable((dimension,dimension),PSD=True)
        solver_map=g if transform is None else g@sp.kron(sp.csr_matrix(transform),sp.csr_matrix(transform),format='csr')
        solver_maps.append(solver_map)
        variables.append(q);maps.append(g);dimensions.append(dimension);transforms.append(transform)
    if variables:
        expression=expression+sp.hstack(solver_maps,format='csr')@cp.hstack([cp.reshape(q,(d*d,),order='C') for q,d in zip(variables,dimensions)])
    objective=b
    if residual_penalty:
        residual=cp.Variable(len(rows))
        weights=np.array([1 if tuple(i for c,i in w if c)==tuple(i for c,i in w if not c) else 2 for w in rows])
        expression=expression+residual
        objective=b-cp.norm1(cp.multiply(weights,residual))
    rhs=np.array([float(h.get(w,0)) for w in rows]);constraint=expression==rhs
    problem=cp.Problem(cp.Maximize(objective),[constraint]);build=time.monotonic()-start;start=time.monotonic()
    problem.solve(solver='CLARABEL',tol_gap_abs=eps,tol_gap_rel=eps,tol_feas=eps,max_iter=200)
    elapsed=time.monotonic()-start
    if b.value is None or x.value is None or constraint.dual_value is None or any(q.value is None for q in variables):
        raise RuntimeError(f'No primal/dual proposal: {problem.status}')
    proposal={'b':float(b.value),'grams':[q.value if tr is None else tr@q.value@tr.T for q,tr in zip(variables,transforms)],'x':np.asarray(x.value),'basis':basis,
              'status':problem.status,'coefficient_rows':len(rows),'gram_dimensions':dimensions,
              'coefficient_map_nonzeros':sum(g.nnz for g in maps),'multiplier_basis_dimension':len(basis),
              'build_seconds':build,'solve_seconds':elapsed}
    proposal['numerical_objective']=float(problem.value)
    factors=[]
    for q,tr in zip(variables,transforms):
        if tr is None:factors.append(None);continue
        eigenvalues,eigenvectors=np.linalg.eigh((q.value+q.value.T)/2)
        factors.append((np.sqrt(np.maximum(eigenvalues,0))[:,None]*eigenvectors.T)@tr.T)
    proposal['factor_matrices']=factors
    return proposal,{'rows':rows,'maps':maps,'multiplier_map':a,'rhs':rhs,'dual':np.asarray(constraint.dual_value),
                     'iterations':problem.solver_stats.num_iters}
