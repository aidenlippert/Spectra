from propose import *
labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');st=time.monotonic()
arr=expm_multiply(-1j*K,psi,start=0,stop=21,num=169,traceA=-1j*K.diagonal().sum());P,s,_=linalg.svd(np.concatenate([psi[:,None],arr.real.T,arr.imag.T],axis=1),full_matrices=False)
np.savez_compressed(ROOT/'development/baseline_basis.npz',J=P,s=s)
for d in [8,12,16,24,32,48]:
 J=P[:,:d];h=J.T@K@J;E=K@J-J@h;z=J.T@psi;ev,V=linalg.eigh(h);ts=np.linspace(0,21,421);tr=(V@(np.exp(-1j*ev[:,None]*ts)*(V.T@z)[:,None])).T
 val=float(np.real(tr[-1].conj()@(J.T@(D[:,None]*J))@tr[-1]));err=np.linalg.norm(psi-J@z)+np.trapezoid(np.sqrt(np.maximum(np.einsum('ti,ij,tj->t',tr.conj(),E.T@E,tr).real,0)),ts)
 print(d,val,err,flush=True)
print('total',time.monotonic()-st)
