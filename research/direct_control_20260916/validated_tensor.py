"""A posteriori tensor error bounds under IEEE binary64 arithmetic.

SVD/QR outputs are untrusted proposals. Their actual local reconstruction
residuals and isometry defects are enclosed. BLAS products use the conservative
standard dot-product bound gamma_(16*k+16) for complex arithmetic, valid for
any parenthesization with at most that many rounding operations per entry.
Finite arithmetic, round-to-nearest, and no flush-to-zero are prerequisites.
All bounds are rounded outwards; subnormal allowances are included.
"""
import math
from fractions import Fraction as F
import numpy as np
from scipy.linalg import svd

U = 2.0**-53
ETA = float.fromhex('0x0.0000000000001p-1022')


def up(x):
    x=float(x)
    if not math.isfinite(x) or x < 0:
        raise ArithmeticError('Nonfinite or negative error bound')
    return math.nextafter(x, math.inf) if x else 0.0


def plus(*xs):
    z=0.0
    for x in xs: z=up(z+x)
    return z


def times(*xs):
    z=1.0
    for x in xs: z=up(z*x)
    return z


def root(x): return up(math.sqrt(up(x)))

def gamma(n):
    x=times(float(n),U)
    if x >= .1: raise ArithmeticError('Rounding budget too large')
    return up(x / math.nextafter(1.0-x, -math.inf))


def frob(a):
    a=np.asarray(a)
    if not np.all(np.isfinite(a)): raise ArithmeticError('Nonfinite tensor')
    # All summands are positive; using a linear reduction bound also covers
    # NumPy's pairwise reduction and a fused implementation.
    s=float(np.sum(a.real*a.real + a.imag*a.imag, dtype=np.float64))
    g=gamma(4*a.size+8)
    return root(plus(up(s / math.nextafter(1-g,-math.inf)), times(4*a.size+8,ETA)))


def mm(a,b):
    c=a@b
    if not np.all(np.isfinite(c)): raise ArithmeticError('Nonfinite product')
    err=plus(times(gamma(16*a.shape[1]+16),frob(a),frob(b)),times(16*a.shape[1]+16,c.size,ETA))
    return c,err


def difference_norm(a,b):
    return plus(frob(a-b),times(gamma(4),plus(frob(a),frob(b))),times(a.size,ETA))


def opnorm(a):
    """Upper bound, not an eigensolver's estimate."""
    a=np.asarray(a)
    if a.shape[0]>=a.shape[1]: g,e=mm(a.conj().T,a)
    else: g,e=mm(a,a.conj().T)
    n=g.shape[0]
    center=float(np.trace(g).real/n)
    shifted=root(plus(abs(center),difference_norm(g,center*np.eye(n)),e))
    ab=np.abs(g)
    rs=up(float(np.max(np.sum(ab,axis=1)))/(1-gamma(4*n+16)))
    cs=up(float(np.max(np.sum(ab,axis=0)))/(1-gamma(4*n+16)))
    induced=root(plus(root(times(rs,cs)),e))
    return min(frob(a),shifted,induced)


def tensor_bounds(tensors):
    left=[];right=[]
    for a in tensors:
        l,s,r=a.shape
        left.append(opnorm(a.reshape(l*s,r)))
        right.append(opnorm(a.reshape(l,s*r)))
    prefix=[1.0]
    for v in left: prefix.append(times(prefix[-1],v))
    suffix=[1.0]*(len(right)+1)
    for i in range(len(right)-1,-1,-1):suffix[i]=times(right[i],suffix[i+1])
    return prefix,suffix


def prefix_gram_bounds(tensors):
    env=np.ones((1,1),complex);err=0.0;bounds=[1.0]
    for a in tensors:
        out=np.zeros((a.shape[2],a.shape[2]),complex);ee=0.0
        for s in range(a.shape[1]):
            p,e1=mm(env,a[:,s,:]);q,e2=mm(a[:,s,:].conj().T,p)
            ee=plus(ee,times(opnorm(a[:,s,:]),e1),e2,
                    times(gamma(4),plus(frob(out),frob(q))))
            out+=q
        on=opnorm(a.reshape(-1,a.shape[2]))
        err=plus(times(on,on,err),ee);env=out
        bounds.append(root(plus(opnorm(env),err)))
    return bounds


def right_canonicalize(tensors,proposal=None):
    """QR proposal with a checked telescoping error; physical dimension arbitrary."""
    pref=prefix_gram_bounds(tensors)
    out=[None]*len(tensors);carry=np.ones((1,1),complex)
    error=0.0;suffix=1.0
    for i in range(len(tensors)-1,-1,-1):
        a=tensors[i];l,p,r=a.shape
        raw,arith=mm(a.reshape(l*p,r),carry)
        raw=raw.reshape(l,p*carry.shape[1])
        if i==0:
            new=raw.reshape(1,p,carry.shape[1]) if proposal is None else proposal[i]
            defect=plus(arith,difference_norm(raw,new.reshape(1,-1)))
        else:
            if proposal is None:
                q,_=np.linalg.qr(raw.T,mode='reduced');c=q.T
            else:c=proposal[i].reshape(proposal[i].shape[0],-1)
            carry=raw@c.conj().T
            rec,re=mm(carry,c)
            defect=plus(arith,difference_norm(raw,rec),re)
            new=c.reshape(c.shape[0],p,-1)
        out[i]=new
        error=plus(error,times(pref[i],defect,suffix))
        suffix=times(suffix,opnorm(new.reshape(new.shape[0],-1)))
    return out,error


def dense_mpo(mpo):
    """Exact rational local entries converted with a charged global HS error."""
    from collections import defaultdict
    tensors=[];local_errors=[]
    for i,layer in enumerate(mpo['layers']):
        exact=defaultdict(F)
        for l,r,mat,c in layer:
            for s,x in enumerate(mat):
                if x:exact[l,s,r]+=F(c)*x
        a=np.zeros((mpo['widths'][i],4,mpo['widths'][i+1]),complex);e2=0.0
        for (l,s,r),c in exact.items():
            v,e=rational_float(c);a[l,s,r]=v;e2=plus(e2,times(e,e))
        tensors.append(a);local_errors.append(root(e2))
    pref,suf=tensor_bounds(tensors)
    exact_suf=[1.0]*(len(tensors)+1)
    for i in range(len(tensors)-1,-1,-1):
        exact_suf[i]=times(plus(opnorm(tensors[i].reshape(tensors[i].shape[0],-1)),local_errors[i]),exact_suf[i+1])
    err=plus(*(times(pref[i],local_errors[i],exact_suf[i+1]) for i in range(len(tensors))))
    return tensors,err


def zip_apply_dense(tensors,operators,max_bond,proposal=None):
    """Canonical dense MPO factors, with no global matrix or product MPS."""
    _,stail=tensor_bounds(tensors);_,otail=tensor_bounds(operators)
    carry=np.ones((1,1,1),complex);output=[];error=0.0;prefix=1.0;local=[]
    for i,(a,w) in enumerate(zip(tensors,operators)):
        chi,wl,dl=carry.shape;dr=a.shape[2];wr=w.shape[2]
        temp,e1=mm(carry.reshape(chi*wl,dl),a.reshape(dl,2*dr))
        temp=temp.reshape(chi,wl,2,dr).transpose(1,2,0,3).reshape(wl*2,chi*dr)
        ww=w.reshape(wl,2,2,wr).transpose(1,3,0,2).reshape(2*wr,wl*2)
        contracted,e2=mm(ww,temp)
        arith=plus(times(opnorm(ww),e1),e2)
        matrix=contracted.reshape(2,wr,chi,dr).transpose(2,0,1,3).reshape(chi*2,wr*dr)
        if i==len(tensors)-1:
            new=matrix.reshape(chi,2,1) if proposal is None else proposal[i]
            defect=plus(arith,difference_norm(matrix,new.reshape(chi*2,1)))
        else:
            if proposal is None:
                u,s,vh=svd(matrix,full_matrices=False,check_finite=False,lapack_driver='gesdd')
                keep=min(max_bond,len(s));u=u[:,:keep]
            else:
                keep=proposal[i].shape[2];u=proposal[i].reshape(chi*2,keep)
                if keep>max_bond:raise ValueError('Proposal bond limit')
            newcarry=u.conj().T@matrix
            rec,re=mm(u,newcarry)
            defect=plus(arith,difference_norm(matrix,rec),re)
            new=u.reshape(chi,2,keep);carry=newcarry.reshape(keep,wr,dr)
        output.append(new)
        # Vectorized operator tail has wl orthogonal-input columns. Its
        # Hilbert-Schmidt norm <= sqrt(wr)*||tail map||; this also bounds
        # the norm of the block operator acting on all MPS right states.
        tail=times(root(wr),otail[i+1],stail[i+1])
        contribution=times(prefix,defect,tail)
        error=plus(error,contribution)
        local.append({'site':i,'local_residual':defect,'tail_upper':tail,'contribution':contribution})
        if i<len(tensors)-1:prefix=times(prefix,opnorm(new.reshape(chi*2,-1)))
    return output,error,local


def rational_float(x):
    x=F(x); y=float(x)
    d=abs(F.from_float(y)-x)
    return y,up(float(d))


def load_mps(cert):
    tensors=[];errors=[]
    for site,edges in enumerate(cert['tensors']):
        a=np.zeros((len(cert['bond_charges'][site]),2,len(cert['bond_charges'][site+1])),complex)
        e2=0.0
        for l,s,r,v in edges:
            x,e=rational_float(F(v,cert['denominator']))
            a[l,s,r]=x;e2=plus(e2,times(e,e))
        tensors.append(a);errors.append(root(e2))
    prefix,suffix=tensor_bounds(tensors)
    # Exact source suffix factors differ from the floating factors by at most
    # the local Frobenius conversion error.
    exact_suffix=[1.0]*(len(tensors)+1)
    for i in range(len(tensors)-1,-1,-1):
        a=tensors[i]
        exact_suffix[i]=times(plus(opnorm(a.reshape(a.shape[0],-1)),errors[i]),exact_suffix[i+1])
    error=plus(*(times(prefix[i],errors[i],exact_suffix[i+1]) for i in range(len(tensors))))
    return tensors,error


def scale_mps(tensors,c):
    out=list(tensors);out[0]=tensors[0]*c
    _,suffix=tensor_bounds(tensors)
    err=times(gamma(8),abs(c),frob(tensors[0]),suffix[1])
    return out,err


def add_mps(a,b,ca=1.0,cb=1.0):
    """Exact block sum except explicitly charged scalar multiplications."""
    aa,ea=scale_mps(a,ca);bb,eb=scale_mps(b,cb)
    out=[]
    for i,(x,y) in enumerate(zip(aa,bb)):
        if i==0:z=np.concatenate((x,y),axis=2)
        elif i==len(a)-1:z=np.concatenate((x,y),axis=0)
        else:
            z=np.zeros((x.shape[0]+y.shape[0],2,x.shape[2]+y.shape[2]),complex)
            z[:x.shape[0],:,:x.shape[2]]=x;z[x.shape[0]:,:,x.shape[2]:]=y
        out.append(z)
    return out,plus(ea,eb)


def mpo_tail_bounds(mpo):
    tails=[None]*(mpo['modes']+1);tails[-1]=[F(1)]
    for i in range(mpo['modes']-1,-1,-1):
        tails[i]=[F(0)]*mpo['widths'][i]
        for l,r,mat,c in mpo['layers'][i]:
            # The constructor's I,Z,creation,annihilation,projector matrices
            # all have spectral norm one. Refuse unsupported matrices.
            if any(type(x)is not int or x not in (-1,0,1) for x in mat):raise ValueError('MPO local factor')
            a,b,c0,d=mat
            # Orthogonal columns or rank-one local factors only.
            col0=a*a+c0*c0;col1=b*b+d*d
            if a*b+c0*d or max(col0,col1)>1:raise ValueError('MPO factor norm > 1')
            tails[i][l]+=abs(F(c))*tails[i+1][r]
    return [root(up(float(sum((x*x for x in t),F())))) for t in tails]


def prepare_mpo(mpo):
    layers=[]
    for layer in mpo['layers']:
        edges=[]
        for l,r,mat,c in layer:
            val,err=rational_float(c)
            for k,x in enumerate(mat):
                if x:edges.append((l,r,k//2,k%2,val*x,err))
        layers.append(edges)
    return layers,mpo_tail_bounds(mpo)


def zip_apply(tensors,mpo,max_bond,prepared=None,proposal=None):
    """Direct MPO action with locally certified compression.

    No MPO*MPS virtual product tensor is built. Only R, one site's contraction,
    and a small-row SVD exist at once. The telescoping certificate charges each
    local residual by a bound on the untouched operator/state tail.
    """
    layers,otail=prepare_mpo(mpo) if prepared is None else prepared
    _,stail=tensor_bounds(tensors)
    carry=np.ones((1,1,1),complex)
    output=[];error=0.0;prefix=1.0;local=[]
    for i,(a,edges) in enumerate(zip(tensors,layers)):
        chi,wl,dl=carry.shape;dr=a.shape[2];wr=mpo['widths'][i+1]
        merged=np.zeros((chi,2,wr,dr),complex)
        caches={};arith=0.0
        for l,r,s,t,c,ce in edges:
            key=(l,t)
            if key not in caches:
                caches[key]=mm(carry[:,l,:],a[:,t,:])
            p,pe=caches[key]
            pn=frob(p)
            term=c*p
            tn=frob(term)
            old=merged[:,s,r,:]
            # Formula coefficients are exact rationals; c is their binary
            # approximation, and ce encloses that conversion error.
            arith=plus(arith,times(abs(c),pe),times(ce,plus(pn,pe)),
                       times(gamma(8),abs(c),pn),
                       times(gamma(4),plus(frob(old),tn)))
            merged[:,s,r,:]=old+term
        matrix=merged.reshape(chi*2,wr*dr)
        if i==len(tensors)-1:
            new=matrix.reshape(chi,2,1) if proposal is None else proposal[i]
            output.append(new)
            defect=plus(arith,difference_norm(matrix,new.reshape(chi*2,1)))
        else:
            if proposal is None:
                u,s,vh=svd(matrix,full_matrices=False,check_finite=False,lapack_driver='gesdd')
                keep=min(max_bond,len(s));u=u[:,:keep]
            else:
                keep=proposal[i].shape[2]
                if keep>max_bond or proposal[i].shape[:2]!=(chi,2):raise ValueError('Proposal bond')
                u=proposal[i].reshape(chi*2,keep)
            # This carry is only an auxiliary witness: its own rounding is
            # immaterial because the actual product U*carry is checked next.
            newcarry=u.conj().T@matrix
            recon,re=mm(u,newcarry)
            defect=plus(arith,difference_norm(matrix,recon),re)
            new=u.reshape(chi,2,keep)
            output.append(new);carry=newcarry.reshape(keep,wr,dr)
        contribution=times(prefix,defect,otail[i+1],stail[i+1])
        error=plus(error,contribution)
        local.append({'site':i,'raw_local_residual_upper':defect,'tail_operator_upper':otail[i+1],
                      'tail_state_upper':stail[i+1],'state_error_contribution':contribution})
        if i<len(tensors)-1:prefix=times(prefix,opnorm(new.reshape(chi*2,-1)))
    return output,error,local


def identity_mpo(m):
    return {'modes':m,'widths':[1]*(m+1),'layers':[[[0,0,[1,0,0,1],'1']] for _ in range(m)]}


def compress(tensors,max_bond):
    return zip_apply(tensors,identity_mpo(len(tensors)),max_bond)


def expectation(tensors,diagonal=None):
    """Expectation center and absolute enclosure radius for a diagonal product.

    Transfer-error propagation uses the positive map's operator norm, bounded
    by ||sum A_s^* A_s||, avoiding a product of all tensor Frobenius norms.
    """
    env=np.ones((1,1),complex);err=0.0
    for i,a in enumerate(tensors):
        weights=(1.,1.) if diagonal is None or i not in diagonal else diagonal[i]
        g=np.zeros((a.shape[2],a.shape[2]),complex);rounderr=0.0
        for s,w in enumerate(weights):
            if not w:continue
            p,e1=mm(env,a[:,s,:]);q,e2=mm(a[:,s,:].conj().T,p)
            e=plus(times(opnorm(a[:,s,:]),e1),e2)
            rounderr=plus(rounderr,times(abs(w),e),times(gamma(8),abs(w),frob(q)),
                          times(gamma(4),plus(frob(g),times(abs(w),frob(q)))))
            g+=w*q
        local_op=opnorm(a.reshape(-1,a.shape[2]))
        gain=times(max(map(abs,weights)),local_op,local_op)
        err=plus(times(gain,err),rounderr);env=g
    if abs(env[0,0].imag)>err:raise ArithmeticError('Hermitian expectation imaginary part')
    return float(env[0,0].real),err
