"""Charge-preserving matrix product states built without a determinant basis."""
from __future__ import annotations
from dataclasses import dataclass
from collections import defaultdict
import numpy as np
from scipy.linalg import svd,qr
PHYS=((0,0),(1,0),(0,1),(1,1))
def plus(a,b):return (a[0]+b[0],a[1]+b[1])
def minus(a,b):return (a[0]-b[0],a[1]-b[1])
def groups(qs):
    out=defaultdict(list)
    for i,q in enumerate(qs):out[tuple(q)].append(i)
    return {q:np.array(ix,dtype=int) for q,ix in out.items()}

@dataclass
class MPS:
    cores:list
    charges:list
    def copy(self):return MPS([a.copy() for a in self.cores],[q.copy() for q in self.charges])
    def scaled(self,value):
        out=self.copy();out.cores[0]=out.cores[0]*value;return out
    @property
    def bond(self):return max(max(a.shape[0],a.shape[2]) for a in self.cores)
    @property
    def entries(self):return sum(a.size for a in self.cores)
    def validate(self):
        if not self.cores or len(self.charges)!=len(self.cores)+1:raise ValueError('chain shape')
        if self.charges[0]!=[(0,0)] or len(self.charges[-1])!=1:raise ValueError('boundaries')
        for i,a in enumerate(self.cores):
            if a.shape!=(len(self.charges[i]),4,len(self.charges[i+1])) or not np.isfinite(a).all():raise ValueError('core shape/value')
            for l,p,r in zip(*np.nonzero(a)):
                if plus(self.charges[i][l],PHYS[p])!=self.charges[i+1][r]:raise ValueError('charge leakage')
        return self

def product(physical):
    cores=[];charges=[[(0,0)]]
    for p in physical:
        if type(p) is not int or not 0<=p<4:raise ValueError('physical state')
        a=np.zeros((1,4,1));a[0,p,0]=1;cores.append(a);charges.append([plus(charges[-1][0],PHYS[p])])
    return MPS(cores,charges)

def overlap(a,b,conjugate=True):
    if len(a.cores)!=len(b.cores):raise ValueError('chain length')
    if a.charges[-1]!=b.charges[-1]:return 0.
    dtype=np.result_type(*a.cores,*b.cores);env={(0,0):np.ones((1,1),dtype=dtype)}
    ga,gb=groups(a.charges[0]),groups(b.charges[0])
    for i,(x,y) in enumerate(zip(a.cores,b.cores)):
        na,nb=groups(a.charges[i+1]),groups(b.charges[i+1]);out={}
        for q,e in env.items():
            for p,ph in enumerate(PHYS):
                qn=plus(q,ph)
                if qn not in na or qn not in nb:continue
                xa=x[:,p,:][np.ix_(ga[q],na[qn])];yb=y[:,p,:][np.ix_(gb[q],nb[qn])]
                if not np.any(xa) or not np.any(yb):continue
                term=(xa.conj().T if conjugate else xa.T)@e@yb
                if qn in out:out[qn]+=term
                else:out[qn]=term
        env=out;ga,gb=na,nb
    return env.get(tuple(a.charges[-1][0]),np.zeros((1,1),dtype=dtype))[0,0]
def norm(a):return float(np.sqrt(max(0.,overlap(a,a).real)))
def amplitude(a,physical):
    if len(physical)!=len(a.cores):raise ValueError('source length')
    v=np.ones(1,dtype=np.result_type(*a.cores))
    for c,p in zip(a.cores,physical):v=v@c[:,p,:]
    return v[0]

def add(a,b,alpha=1.,beta=1.):
    if len(a.cores)!=len(b.cores) or a.charges[-1]!=b.charges[-1]:raise ValueError('sum sector')
    n=len(a.cores);cores=[]
    if n==1:return MPS([alpha*a.cores[0]+beta*b.cores[0]],a.charges)
    for i,(x,y) in enumerate(zip(a.cores,b.cores)):
        if i==0:c=np.concatenate([alpha*x,beta*y],axis=2)
        elif i==n-1:c=np.concatenate([x,y],axis=0)
        else:
            c=np.zeros((x.shape[0]+y.shape[0],4,x.shape[2]+y.shape[2]),dtype=np.result_type(x,y,alpha,beta))
            c[:x.shape[0],:,:x.shape[2]]=x;c[x.shape[0]:,:,x.shape[2]:]=y
        cores.append(c)
    return MPS(cores,[a.charges[0]]+[x+y for x,y in zip(a.charges[1:-1],b.charges[1:-1])]+[a.charges[-1]])

def shift_left(x,i):
    a=x.cores[i];l,_,r=a.shape;rg=groups(x.charges[i]);cg=groups([minus(q,ph) for ph in PHYS for q in x.charges[i+1]])
    blocks=[];qs=[]
    for q,rows in rg.items():
        cols=cg.get(q)
        if cols is None:continue
        mat=a.reshape(l,4*r)[np.ix_(rows,cols)]
        if not np.any(mat):continue
        Q,R=qr(mat.T,mode='economic',check_finite=False);blocks.append((rows,cols,Q.T,R.T));qs.extend([q]*Q.shape[1])
    if not qs:raise ArithmeticError('zero tensor in canonicalization')
    core=np.zeros((len(qs),4*r),dtype=a.dtype);transfer=np.zeros((l,len(qs)),dtype=a.dtype);j=0
    for rows,cols,Q,R in blocks:
        k=Q.shape[0];ix=np.arange(j,j+k);core[np.ix_(ix,cols)]=Q;transfer[np.ix_(rows,ix)]=R;j+=k
    b=x.cores[i-1];x.cores[i]=core.reshape(len(qs),4,r);x.charges[i]=qs
    x.cores[i-1]=(b.reshape(-1,l)@transfer).reshape(b.shape[0],4,len(qs))

def shift_right(x,i):
    a=x.cores[i];l,_,r=a.shape;rg=groups([plus(q,ph) for q in x.charges[i] for ph in PHYS]);cg=groups(x.charges[i+1]);blocks=[];qs=[]
    for q,cols in cg.items():
        rows=rg.get(q)
        if rows is None:continue
        mat=a.reshape(l*4,r)[np.ix_(rows,cols)]
        if not np.any(mat):continue
        Q,R=qr(mat,mode='economic',check_finite=False);blocks.append((rows,cols,Q,R));qs.extend([q]*Q.shape[1])
    if not qs:raise ArithmeticError('zero tensor in canonicalization')
    core=np.zeros((l*4,len(qs)),dtype=a.dtype);transfer=np.zeros((len(qs),r),dtype=a.dtype);j=0
    for rows,cols,Q,R in blocks:
        k=Q.shape[1];ix=np.arange(j,j+k);core[np.ix_(rows,ix)]=Q;transfer[np.ix_(ix,cols)]=R;j+=k
    b=x.cores[i+1];x.cores[i]=core.reshape(l,4,len(qs));x.charges[i+1]=qs
    x.cores[i+1]=(transfer@b.reshape(r,-1)).reshape(len(qs),4,b.shape[2])

def right_canonical(a):
    x=a.copy()
    for i in range(len(x.cores)-1,0,-1):shift_left(x,i)
    return x

def split(mat,rowq,colq,max_bond,cutoff):
    rg,cg=groups(rowq),groups(colq);blocks=[];sing=[]
    for q,cols in cg.items():
        rows=rg.get(q)
        if rows is None:continue
        part=mat[np.ix_(rows,cols)]
        if not np.any(part):continue
        U,S,Vh=svd(part,full_matrices=False,check_finite=False,lapack_driver='gesdd');idx=len(blocks)
        blocks.append((q,rows,cols,U,S,Vh));sing.extend((float(v),idx,j) for j,v in enumerate(S))
    sing.sort(reverse=True);selected=sing[:max_bond]
    while len(selected)>1 and selected[-1][0]<=cutoff:selected.pop()
    if not selected or not selected[0][0]:raise ArithmeticError('zero split')
    keep=defaultdict(list)
    for _,idx,j in selected:keep[idx].append(j)
    rank=len(selected);core=np.zeros((mat.shape[0],rank),dtype=mat.dtype);transfer=np.zeros((rank,mat.shape[1]),dtype=mat.dtype);qs=[];j=0
    for idx,(q,rows,cols,U,S,Vh) in enumerate(blocks):
        js=sorted(keep.get(idx,[]));k=len(js)
        if not k:continue
        ix=np.arange(j,j+k);core[np.ix_(rows,ix)]=U[:,js];transfer[np.ix_(ix,cols)]=S[js,None]*Vh[js,:];qs.extend([q]*k);j+=k
    return core,transfer,qs,sum(v*v for v,_,_ in sing[len(selected):])

def compress(a,max_bond=64,cutoff=1e-11):
    if type(max_bond) is not int or max_bond<1 or cutoff<0:raise ValueError('compression parameters')
    out=right_canonical(a);discard=0.
    for i in range(len(out.cores)-1):
        x=out.cores[i];l,_,r=x.shape
        core,transfer,qs,d=split(x.reshape(l*4,r),[plus(q,ph) for q in out.charges[i] for ph in PHYS],out.charges[i+1],max_bond,cutoff)
        out.cores[i]=core.reshape(l,4,len(qs));nxt=out.cores[i+1];out.cores[i+1]=(transfer@nxt.reshape(r,-1)).reshape(len(qs),4,nxt.shape[2]);out.charges[i+1]=qs;discard+=d
    return out,float(np.sqrt(discard))

def apply(mpo,ket):
    """Uncompressed local action, reserved for small tests."""
    out=[];charges=[[plus(q,d) for d in ds for q in qs] for ds,qs in zip(mpo['charges'],ket.charges)]
    for ops,a,left,right in zip(mpo['cores'],ket.cores,mpo['charges'][:-1],mpo['charges'][1:]):
        dl,_,dr=a.shape;c=np.zeros((len(left)*dl,4,len(right)*dr),dtype=np.result_type(a,float))
        for (u,v),op in ops.items():
            for p,q in zip(*np.nonzero(op)):c[u*dl:(u+1)*dl,p,v*dr:(v+1)*dr]+=op[p,q]*a[:,q,:]
        out.append(c)
    return MPS(out,charges)

def apply_compressed(mpo,ket,max_bond=64,cutoff=1e-11):
    """Streaming operator action with local charge SVDs. Discard is diagnostic."""
    cores=[];charges=[[(0,0)]];carry=np.ones((1,1),dtype=np.result_type(*ket.cores));discard=0.;largest=0
    for i,(ops,a) in enumerate(zip(mpo['cores'],ket.cores)):
        dl,_,dr=a.shape;wr=len(mpo['charges'][i+1]);nl=carry.shape[0];c=np.zeros((nl,4,wr*dr),dtype=np.result_type(a,carry,float))
        for (u,v),op in ops.items():
            left=carry[:,u*dl:(u+1)*dl]
            for p,q in zip(*np.nonzero(op)):c[:,p,v*dr:(v+1)*dr]+=op[p,q]*(left@a[:,q,:])
        largest=max(largest,c.size);rq=[plus(q,d) for d in mpo['charges'][i+1] for q in ket.charges[i+1]]
        if i==len(ket.cores)-1:cores.append(c);charges.append(rq);break
        try:core,carry,qs,d=split(c.reshape(nl*4,wr*dr),[plus(q,ph) for q in charges[-1] for ph in PHYS],rq,max_bond,cutoff)
        except ArithmeticError:return ket.scaled(0.),dict(max_temporary_entries=largest,discard_diagnostic=float(np.sqrt(discard)))
        cores.append(core.reshape(nl,4,len(qs)));charges.append(qs);discard+=d
    return MPS(cores,charges),dict(max_temporary_entries=largest,discard_diagnostic=float(np.sqrt(discard)))
