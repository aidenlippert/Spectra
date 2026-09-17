"""Numerical block-QR/block-SVD MPS with exactly enforced charge patterns."""
import numpy as np
from research.constructive_compression_20260916 import mps_numeric as mn

def plus(q,r):return (q[0]+r[0],q[1]+r[1])
def occ(i,s):return (s,0) if i%2==0 else (0,s)
def _check(a,q):
    mn._check(a)
    if len(q)!=len(a)+1:raise ValueError('Charge length')
    if any(x.shape!=(len(q[i]),2,len(q[i+1])) for i,x in enumerate(a)):raise ValueError('Charge shape')

def structural_error(arrays,charges):
    worst=0.
    for i,a in enumerate(arrays):
        for l,ql in enumerate(charges[i]):
            for s in (0,1):
                for r,qr in enumerate(charges[i+1]):
                    if plus(ql,occ(i,s))!=tuple(qr):worst=max(worst,abs(a[l,s,r]))
    return worst

def linear_combination(states,coeffs):
    for a,q in states:_check(a,q)
    if any(q[0]!=states[0][1][0] or q[-1]!=states[0][1][-1] for a,q in states):raise ValueError('Different sectors')
    arrays=mn.linear_combination([a for a,q in states],coeffs)
    L=len(arrays);q=[list(states[0][1][0])]
    for i in range(1,L):q.append([tuple(x) for a,qq in states for x in qq[i]])
    q.append(list(states[0][1][-1]));return arrays,q

def apply_mpo(arrays,charges,W,operator_charges):
    _check(arrays,charges)
    out=mn.apply_mpo(arrays,W)
    qs=[[plus(qo,qs) for qo in operator_charges[i] for qs in charges[i]] for i in range(len(charges))]
    return out,qs

def compress(arrays,charges,max_bond):
    _check(arrays,charges)
    if type(max_bond)!=int or max_bond<1:raise ValueError('Positive bond cap')
    out=[a.copy() for a in arrays];qs=[[tuple(x) for x in q] for q in charges];discarded=0.
    # Right-canonicalize independently in every charge block.
    for i in range(len(out)-1,0,-1):
        l,p,r=out[i].shape;raw=out[i].reshape(l,p*r);pieces=[]
        for q in sorted(set(qs[i])):
            rr=[j for j,x in enumerate(qs[i]) if x==q]
            cc=[s*r+j for s in (0,1) for j,x in enumerate(qs[i+1]) if plus(q,occ(i,s))==x]
            if not cc:continue
            qmat,rmat=np.linalg.qr(raw[np.ix_(rr,cc)].T,mode='reduced')
            pieces.append((q,rr,cc,qmat.T,rmat.T))
        rank=sum(x[3].shape[0] for x in pieces)
        if not rank:raise ValueError('No charge-compatible QR block')
        new=np.zeros((rank,p*r));carry=np.zeros((l,rank));labels=[];at=0
        for q,rr,cc,b,c in pieces:
            k=b.shape[0];new[np.ix_(range(at,at+k),cc)]=b;carry[np.ix_(rr,range(at,at+k))]=c
            labels.extend([q]*k);at+=k
        out[i]=new.reshape(rank,p,r);out[i-1]=np.einsum('aps,sr->apr',out[i-1],carry,optimize=True);qs[i]=labels
    # TT-SVD with one total bond cap across all charge blocks.
    for i in range(len(out)-1):
        l,p,r=out[i].shape;raw=out[i].reshape(l*p,r);pieces={};candidates=[]
        for q in sorted(set(qs[i+1])):
            rr=[a*p+s for a,ql in enumerate(qs[i]) for s in (0,1) if plus(ql,occ(i,s))==q]
            cc=[b for b,qr in enumerate(qs[i+1]) if qr==q]
            if not rr or not cc:continue
            u,s,v=np.linalg.svd(raw[np.ix_(rr,cc)],full_matrices=False);pieces[q]=(rr,cc,u,s,v)
            candidates.extend((float(value),q,k) for k,value in enumerate(s) if value>1e-14)
            discarded+=float(np.sum(s[s<=1e-14]**2))
        candidates.sort(key=lambda z:z[0],reverse=True);chosen=candidates[:max_bond]
        discarded+=sum(x[0]**2 for x in candidates[max_bond:]);chosen.sort(key=lambda z:(z[1],z[2]))
        if not chosen:raise ValueError('Numerical zero state')
        new=np.zeros((l*p,len(chosen)));carry=np.zeros((len(chosen),r));labels=[]
        for at,(value,q,k) in enumerate(chosen):
            rr,cc,u,s,v=pieces[q];new[rr,at]=u[:,k];carry[at,cc]=s[k]*v[k];labels.append(q)
        out[i]=new.reshape(l,p,len(chosen));out[i+1]=np.einsum('ab,bpc->apc',carry,out[i+1],optimize=True);qs[i+1]=labels
    return out,qs,{'discarded_squared_norm':discarded,'bond_dimensions':[len(x) for x in qs]}
