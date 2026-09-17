from propose import *
started=time.monotonic();labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');Dp=sparse.diags(D)
sn=[psi[:,None]]
for u,v in [(.2,.1),(.202,.098),(.198,.102)]:
 A=K+u*Dp+v*W;arr=expm_multiply(-1j*A,psi,start=0,stop=21,num=169,traceA=-1j*A.diagonal().sum());sn.extend([arr.real.T,arr.imag.T])
X=np.concatenate(sn,axis=1);P,s,_=linalg.svd(X,full_matrices=False,check_finite=False)
np.savez_compressed(ROOT/'development/long_adapted_basis.npz',J=P,s=s)
print('spectrum',[(r,float(s[r])) for r in [8,16,24,32,48,64,96,128]],flush=True)
for d in [8,12,16,24,32,48,64,96]:
 J=P[:,:d];hs=np.array([J.T@K@J,J.T@(D[:,None]*J),J.T@W@J]);E0=K@J-J@hs[0];ED=D[:,None]*J-J@hs[1];EW=W@J-J@hs[2]
 z=J.T@psi;err=np.linalg.norm(psi-J@z);A=hs[0]+.2*hs[1]+.1*hs[2];ev,V=linalg.eigh(A);grid=np.linspace(0,21,421)
 arr=(V@(np.exp(-1j*ev[:,None]*grid)*(V.T@z)[:,None])).T;z=arr[-1];E=E0+.2*ED+.1*EW;R=E.T@E
 vals=np.maximum(np.einsum('ti,ij,tj->t',arr.conj(),R,arr).real,0);err+=np.trapezoid(np.sqrt(vals),grid)
 print('d',d,'D',float(np.real(z.conj()@hs[1]@z)),'err_estimate',err,flush=True)
print('total',time.monotonic()-started)
