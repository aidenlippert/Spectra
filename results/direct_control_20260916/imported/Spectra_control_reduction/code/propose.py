import os
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'
import json, time, itertools, sys
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy import sparse, linalg
from scipy.sparse.linalg import expm_multiply
from numba import njit
ROOT=Path(__file__).resolve().parents[1]
@njit
def build(labels, words, cs, lens):
    idx=np.full(65536,-1,np.int32)
    for i in range(len(labels)):idx[labels[i]]=i
    rows=np.empty(len(labels)*len(words),np.int32);cols=np.empty_like(rows);data=np.empty(len(rows),np.int64);nn=0
    for j in range(len(labels)):
      s0=labels[j]
      for k in range(len(words)):
        s=s0;sg=1
        for p in range(lens[k]-1,-1,-1):
          c,m=words[k,p,0],words[k,p,1]
          if ((s>>m)&1)==c:sg=0;break
          z=s&((1<<m)-1);cnt=0
          while z:cnt+=1;z&=z-1
          if cnt%2:sg=-sg
          s^=(1<<m)
        if sg:
          ii=idx[s]
          if ii<0: raise ValueError('outside fixed spin')
          rows[nn]=ii;cols[nn]=j;data[nn]=sg*cs[k];nn+=1
    return rows[:nn],cols[:nn],data[:nn]
def numeric_mps(state,labels):
    tens=[]
    for i,tt in enumerate(state['tensors']):
      a=np.zeros((len(state['bond_charges'][i]),2,len(state['bond_charges'][i+1])))
      for l,s,r,v in tt:a[l,s,r]=v/state['denominator']
      tens.append(a)
    out=[]
    for lab in labels:
      v=np.array([1.])
      for i,tt in enumerate(tens):v=v@tt[:,(lab>>i)&1,:]
      out.append(v[0])
    return np.array(out)
def setup():
    p=ROOT/'development'/'numeric_inputs.npz'
    if p.exists():
      ar=np.load(p);H=sparse.load_npz(ROOT/'development'/'H.npz');W=sparse.load_npz(ROOT/'development'/'W.npz')
      return ar['labels'],H,ar['D'],W,ar['psi']
    data=json.loads((ROOT/'inputs/fixture.json').read_text());st=json.loads((ROOT/'inputs/state.json').read_text())
    labels=sorted(sum(1<<(2*i) for i in a)+sum(1<<(2*j+1) for j in b) for a in itertools.combinations(range(8),4) for b in itertools.combinations(range(8),4))
    labels=np.array(labels,dtype=np.int64)
    wl=data['hamiltonian'];words=np.zeros((len(wl),4,2),np.int32);lens=np.array([len(e['word']) for e in wl],np.int32)
    cs=np.array([int(F(e['coefficient'])*10**12) for e in wl],np.int64)
    for k,e in enumerate(wl): words[k,:len(e['word']),:]=e['word']
    r,c,v=build(labels,words,cs,lens);HI=sparse.coo_matrix((v,(r,c)),shape=(len(labels),len(labels))).tocsr(); HI.eliminate_zeros()
    assert (HI-HI.T).nnz==0
    sparse.save_npz(ROOT/'development'/'H_integer.npz',HI)
    H=HI.astype(float)/1e12
    D=np.array([((q>>6)&1)+((q>>7)&1)-((q>>10)&1)-((q>>11)&1) for q in labels],float)
    ww=np.array([[[1,i],[0,j],[0,0],[0,0]] for i,j in [(6,10),(10,6),(7,11),(11,7)]],np.int32)
    r,c,v=build(labels,ww,np.ones(4,np.int64),np.full(4,2,np.int32));W=sparse.coo_matrix((v,(r,c)),shape=H.shape).tocsr()
    psi=numeric_mps(st,labels);psi/=np.linalg.norm(psi)
    sparse.save_npz(ROOT/'development'/'H.npz',H);sparse.save_npz(ROOT/'development'/'W.npz',W)
    np.savez_compressed(p,labels=labels,D=D,psi=psi)
    return labels,H,D,W,psi

def discover():
    started=time.monotonic();labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');Dp=sparse.diags(D)
    print('Setup',len(labels),H.nnz,'E',psi@H@psi,'D',psi@(D*psi),time.monotonic()-started,flush=True)
    snapshots=[psi[:,None]]
    training=[(-.1,-.1),(-.1,.1),(.1,-.1),(.1,.1),(0,.15),(0,-.15)]
    for u,v in training:
      t=time.monotonic(); A=K+u*Dp+v*W
      arr=expm_multiply(-1j*A,psi,start=0,stop=20,num=41,traceA=-1j*A.diagonal().sum())
      snapshots.extend([arr.real.T,arr.imag.T]);print('snapshot',u,v,time.monotonic()-t,flush=True)
    X=np.concatenate(snapshots,axis=1);P,s,_=linalg.svd(X,full_matrices=False,check_finite=False)
    print('Svals',[(r,float(s[r])) for r in [8,16,24,32,48,64,96,128,160] if r<len(s)],flush=True)
    np.savez_compressed(ROOT/'development'/'basis_proposal.npz',J=P,s=s)
    # test a new switching history, not in the snapshot training
    controls=[(.07,.12),(-.05,.08),(.03,-.11),(-.08,.15)]
    for d in [16,24,32,48,64,96,128]:
      J=P[:,:d];h=[J.T@K@J,J.T@(D[:,None]*J),J.T@W@J]
      Es=[K@J-J@h[0],D[:,None]*J-J@h[1],W@J-J@h[2]]
      z=J.T@psi;pp=psi.astype(complex);eta=np.linalg.norm(psi-J@z);history=[]
      for u,v in controls:
        hu=h[0]+u*h[1]+v*h[2];ev,V=linalg.eigh(hu)
        grid=np.linspace(0,5,101); traj=(V@(np.exp(-1j*ev[:,None]*grid)*(V.T@z)[:,None])).T
        E=Es[0]+u*Es[1]+v*Es[2];Dg=E.T@E
        res=np.sqrt(np.maximum(np.einsum('ti,ij,tj->t',traj.conj(),Dg,traj).real,0))
        eta+=np.trapezoid(res,grid)
        z=traj[-1];AA=K+u*Dp+v*W
        pp=expm_multiply(-1j*5*AA,pp,traceA=-1j*5*AA.diagonal().sum())
      print('d',d,'a posteriori estimate',eta,'actual err',np.linalg.norm(pp-J@z),'D',float((J@z).conj()@(D*(J@z))).real,flush=True)
    print('total',time.monotonic()-started,flush=True)
if __name__=='__main__':discover()
