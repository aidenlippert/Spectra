"""Numerical polynomial Gram maps built from an exact monomial CAR map.

Only the discovery matrix is floating point. Polynomial generators and the
final SOS exporter remain exact; this module cannot certify an energy bound.
All support word pairs and dense contraction entries are explicitly counted.
"""
import numpy as np
from scipy import sparse
from experiments.marginal_symbolic import word_product


def prepare(groups):
    blocks=[];allwords=set();pairs=0;nnz=0
    for group in groups:
        words=sorted({w for p in group['polynomials'] for w in p},key=lambda w:(len(w),w))
        terms=[];size=len(words);pairs+=size*size
        for i,left in enumerate(words):
            dagger=tuple((1-c,j) for c,j in reversed(left))
            for j,right in enumerate(words):
                for w,c in word_product(dagger,right):
                    terms.append((w,i*size+j,c));allwords.add(w)
        nnz+=len(terms);blocks.append({'words':words,'terms':terms})
    return blocks,allwords,{'monomial_word_pairs':pairs,'monomial_map_nonzeros':nnz}


def contract(block,polynomials,lookup,batch=16,transform=None):
    if type(batch) is not int or batch<1:raise ValueError('Positive contraction batch required')
    words=block['words'];s=len(words);k=len(polynomials)
    terms=block['terms'];g=sparse.csr_matrix(([float(c) for _,_,c in terms],
        ([lookup[w] for w,_,_ in terms],[j for _,j,_ in terms])),shape=(len(lookup),s*s))
    coefficients=np.array([[float(p.get(w,0)) for p in polynomials] for w in words])
    transform_work=0
    if transform is not None:
        transform=np.asarray(transform,dtype=float)
        if transform.shape!=(k,k) or not np.all(np.isfinite(transform)):
            raise ValueError('Invalid polynomial basis transform')
        coefficients=coefficients@transform.T;transform_work=s*k*k
    indices=[(i,j) for i in range(k) for j in range(i,k)]
    chunks=[];peak=0
    for offset in range(0,len(indices),batch):
        pairs=indices[offset:offset+batch];ii=np.array([i for i,j in pairs]);jj=np.array([j for i,j in pairs])
        left=coefficients[:,ii];right=coefficients[:,jj]
        values=left[:,None,:]*right[None,:,:]
        values+=right[:,None,:]*left[None,:,:]*(ii!=jj)[None,None,:]
        transformed=g@values.reshape(s*s,len(pairs))
        if not np.all(np.isfinite(transformed)):raise ValueError('Nonfinite discovery map')
        peak=max(peak,values.nbytes+transformed.nbytes)
        chunks.append(sparse.csc_matrix(transformed))
    matrix=sparse.hstack(chunks,format='csc')
    return matrix,[i*k+j for i,j in indices],{'basis_transform_multiply_entries':transform_work,'dense_outer_entries':2*s*s*len(indices),
        'dense_sum_entries':s*s*len(indices),'sparse_dense_multiply_entries':g.nnz*len(indices),
        'contraction_peak_primary_array_bytes':peak,'contraction_batches':len(chunks),
        'scope':'Floating proposal map only; exact final polynomial replay is mandatory.'}


def whitening_transform(polynomials):
    """Floating square basis change; no directions dropped or bound certified."""
    import time
    start=time.monotonic()
    words=sorted({w for p in polynomials for w in p},key=lambda w:(len(w),w))
    C=np.array([[float(p.get(w,0)) for p in polynomials] for w in words])
    if C.shape[0]<C.shape[1] or not np.all(np.isfinite(C)):
        raise ValueError('Cannot whiten an overcomplete/nonfinite polynomial basis')
    # Disjoint exact monomial supports are orthogonal components. Keep their
    # off-component zeros exact instead of allowing SVD roundoff to densify maps.
    k=len(polynomials);parent=list(range(k))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    owners={}
    for i,p in enumerate(polynomials):
        for w,c in p.items():
            if not c:continue
            if w in owners:parent[find(i)]=find(owners[w])
            else:owners[w]=i
    components={}
    for i in range(k):components.setdefault(find(i),[]).append(i)
    W=np.zeros((k,k));singular=[]
    for ix in components.values():
        active=np.flatnonzero(np.any(C[:,ix]!=0,axis=1))
        block=C[np.ix_(active,ix)]
        _,sv,vt=np.linalg.svd(block,full_matrices=False)
        if len(sv)!=len(ix) or sv[-1]<=sv[0]*1e-12:
            raise ValueError('Numerically rank deficient basis; no automatic truncation')
        W[np.ix_(ix,ix)]=vt/sv[:,None];singular.extend(sv)
    transformed=C@W.T
    err=float(np.max(np.abs(transformed.T@transformed-np.eye(k))))
    if not np.all(np.isfinite(W)) or err>1e-7:
        raise ValueError('Failed numerical whitening check')
    return W,{'basis_dimension':C.shape[1],'support_monomials':C.shape[0],
              'support_component_dimensions':[len(ix) for ix in components.values()],
              'coefficient_matrix_entries':int(C.size),'estimated_condition_before':float(max(singular)/min(singular)),
              'orthogonality_max_error':err,'max_transform_entry':float(np.max(np.abs(W))),
              'seconds':time.monotonic()-start,'scope':'Numerical congruence only; original exact generators retained for export.'}
