"""Numerical real MPS proposals; no global vectors and no acceptance claims."""
import numpy as np

def _check(state):
    if not state or any(x.ndim!=3 for x in state):raise ValueError('MPS tensors')
    if state[0].shape[0]!=1 or state[-1].shape[2]!=1:raise ValueError('Scalar boundaries')
    if any(a.shape[2]!=b.shape[0] for a,b in zip(state,state[1:])):raise ValueError('Bond mismatch')

def inner(a,b):
    _check(a);_check(b)
    if len(a)!=len(b):raise ValueError('Length mismatch')
    e=np.ones((1,1))
    for x,y in zip(a,b):
        if x.shape[1]!=y.shape[1]:raise ValueError('Physical dimensions')
        e=np.einsum('ab,asi,bsj->ij',e,x.conj(),y,optimize=True)
    return float(e[0,0].real)

def norm(a):return float(np.sqrt(max(0.,inner(a,a))))

def linear_combination(states,coeffs):
    if not states or len(states)!=len(coeffs):raise ValueError('Empty or misaligned combination')
    for x in states:_check(x)
    L=len(states[0]);out=[]
    if any(len(x)!=L for x in states):raise ValueError('Length mismatch')
    if L==1:return [sum(c*x[0] for c,x in zip(coeffs,states))]
    for i in range(L):
        p=states[0][i].shape[1]
        if any(x[i].shape[1]!=p for x in states):raise ValueError('Physical dimensions')
        if i==0:t=np.concatenate([c*x[i] for c,x in zip(coeffs,states)],axis=2)
        elif i==L-1:t=np.concatenate([x[i] for x in states],axis=0)
        else:
            t=np.zeros((sum(x[i].shape[0] for x in states),p,sum(x[i].shape[2] for x in states)))
            a=b=0
            for x in states:
                l,_,r=x[i].shape;t[a:a+l,:,b:b+r]=x[i];a+=l;b+=r
        out.append(t)
    return out

def apply_mpo(state,W):
    _check(state)
    if len(state)!=len(W):raise ValueError('Length mismatch')
    out=[]
    for a,w in zip(state,W):
        if w.ndim!=4 or a.shape[1]!=w.shape[3]:raise ValueError('MPO shape')
        t=np.einsum('abpk,ikj->aipbj',w,a,optimize=True)
        out.append(t.reshape(w.shape[0]*a.shape[0],w.shape[2],w.shape[1]*a.shape[2]))
    return out

def compress(state,max_bond):
    _check(state)
    if type(max_bond)!=int or max_bond<1:raise ValueError('Positive bond cap')
    out=[a.copy() for a in state]
    for i in range(len(out)-1,0,-1):
        l,p,r=out[i].shape;q,rr=np.linalg.qr(out[i].reshape(l,p*r).T,mode='reduced')
        out[i]=q.T.reshape(q.shape[1],p,r)
        out[i-1]=np.einsum('aps,sr->apr',out[i-1],rr.T,optimize=True)
    discarded=0.
    for i in range(len(out)-1):
        l,p,r=out[i].shape;u,s,v=np.linalg.svd(out[i].reshape(l*p,r),full_matrices=False)
        k=min(max_bond,len(s));discarded+=float(s[k:]@s[k:])
        out[i]=u[:,:k].reshape(l,p,k)
        out[i+1]=np.einsum('ab,bpc->apc',s[:k,None]*v[:k],out[i+1],optimize=True)
    return out,{'discarded_squared_norm':discarded,'bond_dimensions':[a.shape[2] for a in out]}
