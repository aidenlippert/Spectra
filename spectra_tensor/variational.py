"""Local response equations in adaptive charge-conserving tensor spaces.

Completed sweeps are not success proofs. The exact final residual determines
acceptance even when some local iterative equations did not fully converge.
"""
from time import perf_counter
import numpy as np
from scipy.sparse.linalg import LinearOperator,gmres
from . import tensor as t,physics as p
from .solver import screen

def left_update(env,a,ops,size):
    out=[np.zeros((a.shape[2],a.shape[2]),complex) for _ in range(size)]
    for (u,v),op in ops.items():
        for pp,q in zip(*np.nonzero(op)):out[v]+=op[pp,q]*(a[:,pp,:].conj().T@env[u]@a[:,q,:])
    return out

def right_update(env,a,ops,size):
    out=[np.zeros((a.shape[0],a.shape[0]),complex) for _ in range(size)]
    for (u,v),op in ops.items():
        for pp,q in zip(*np.nonzero(op)):out[u]+=op[pp,q]*(a[:,pp,:].conj()@env[v]@a[:,q,:].T)
    return out

def local_solve(x,i,L,R,bl,br,ops,dl,dr,z,source,rtol=1e-6,maxiter=4):
    gl,gr=t.groups(x.charges[i]),t.groups(x.charges[i+1]);layout={};size=0
    for q,rows in gl.items():
        for pp,ph in enumerate(t.PHYS):
            qn=t.plus(q,ph)
            if qn not in gr:continue
            cols=gr[qn];count=len(rows)*len(cols);layout[q,pp]=(slice(size,size+count),(len(rows),len(cols)),rows,cols);size+=count
    if not size:raise ArithmeticError('empty local space')
    def unpack(v):return {key:v[sl].reshape(shape) for key,(sl,shape,_,_) in layout.items()}
    terms=[]
    for (u,v),op in ops.items():
        for pp,qq in zip(*np.nonzero(op)):
            for ql in gl:
                ink=(ql,qq);outk=(t.plus(ql,dl[u]),pp)
                if ink not in layout or outk not in layout:continue
                if t.plus(t.plus(ql,t.PHYS[qq]),dr[v])!=t.plus(outk[0],t.PHYS[pp]):continue
                _,_,ci,di=layout[ink];_,_,ai,bi=layout[outk];lb=L[u][np.ix_(ai,ci)];rb=R[v][np.ix_(bi,di)].T
                if np.any(lb) and np.any(rb):terms.append((ink,outk,op[pp,qq]*lb,rb))
    def h(v):
        inp=unpack(v);out=np.zeros(size,complex);dest=unpack(out)
        for ink,outk,l,r in terms:dest[outk]+=l@inp[ink]@r
        return out
    calls=[0]
    def action(v):calls[0]+=1;return z*v-h(v)
    rhs=np.zeros(size,complex);guess=rhs.copy();diag=rhs.copy();db=unpack(diag)
    for (q,pp),(sl,shape,rows,cols) in layout.items():
        guess[sl]=x.cores[i][:,pp,:][np.ix_(rows,cols)].reshape(-1)
        if pp==source:rhs[sl]=(bl[rows,None]*br[None,cols]).reshape(-1)
    for ink,outk,l,r in terms:
        if ink==outk:db[ink]+=np.diag(l)[:,None]*np.diag(r)[None,:]
    pre=z-diag;A=LinearOperator((size,size),matvec=action,dtype=complex);M=LinearOperator((size,size),matvec=lambda v:v/pre,dtype=complex)
    sol,info=gmres(A,rhs,x0=guess,M=M,rtol=rtol,atol=1e-12,maxiter=maxiter,restart=min(32,size))
    a=np.zeros_like(x.cores[i],dtype=complex)
    for (q,pp),(sl,shape,rows,cols) in layout.items():a[:,pp,:][np.ix_(rows,cols)]=sol[sl].reshape(shape)
    x.cores[i]=a;return dict(local_variables=size,matvecs=calls[0],info=int(info))

def sweep_solve(spec,z,initial,sweeps=4,target=.00055,verbose=True):
    if sweeps<1 or target<=0:raise ValueError('sweep parameters')
    start=perf_counter();mpo=p.compile_hubbard(spec);source=p.source(spec);n=len(source);x=t.right_canonical(initial);history=[];total=0;best=None
    for sweep in range(sweeps):
        R=[None]*(n+1);br=[None]*(n+1);R[n]=[np.ones((1,1),complex)];br[n]=np.ones(1,complex)
        for i in range(n-1,-1,-1):
            R[i]=right_update(R[i+1],x.cores[i],mpo['cores'][i],len(mpo['charges'][i]));br[i]=x.cores[i][:,source[i],:].conj()@br[i+1]
        L=[None]*(n+1);bl=[None]*(n+1);L[0]=[np.ones((1,1),complex)];bl[0]=np.ones(1,complex);largest=failed=0
        for i in range(n):
            stats=local_solve(x,i,L[i],R[i+1],bl[i],br[i+1],mpo['cores'][i],mpo['charges'][i],mpo['charges'][i+1],z,source[i]);total+=stats['matvecs'];failed+=stats['info']!=0;largest=max(largest,stats['local_variables'])
            if i<n-1:t.shift_right(x,i)
            L[i+1]=left_update(L[i],x.cores[i],mpo['cores'][i],len(mpo['charges'][i+1]));bl[i+1]=x.cores[i][:,source[i],:].conj().T@bl[i]
        R[n]=[np.ones((1,1),complex)];br[n]=np.ones(1,complex)
        for i in range(n-1,-1,-1):
            stats=local_solve(x,i,L[i],R[i+1],bl[i],br[i+1],mpo['cores'][i],mpo['charges'][i],mpo['charges'][i+1],z,source[i]);total+=stats['matvecs'];failed+=stats['info']!=0;largest=max(largest,stats['local_variables'])
            if i>0:t.shift_left(x,i)
            R[i]=right_update(R[i+1],x.cores[i],mpo['cores'][i],len(mpo['charges'][i]));br[i]=x.cores[i][:,source[i],:].conj()@br[i+1]
        result=screen(spec,x,z);row=dict(sweep=sweep+1,radius=result['radius'],seconds=perf_counter()-start,bond=x.bond,max_local_variables=largest,local_matvecs=total,local_nonconvergence=failed);history.append(row)
        if verbose:print(row,flush=True)
        if best is None or result['radius']<best[0]:best=(result['radius'],x.copy())
        if result['radius']<=target:break
    return best[1],dict(history=history,seconds=perf_counter()-start,enumerated_states=0,local_matvecs=total)

def enrich(spec,z,x,max_bond,amplitude=.01):
    hx,stats=t.apply_compressed(p.compile_hubbard(spec),x,max_bond,1e-12)
    trial=t.add(x,hx,alpha=1-amplitude*z,beta=amplitude);trial=t.add(trial,t.product(p.source(spec)),beta=amplitude);trial,discard=t.compress(trial,max_bond,1e-12)
    trial.validate();return trial,dict(operator_action=stats,discard_diagnostic=discard,new_bond=trial.bond)
