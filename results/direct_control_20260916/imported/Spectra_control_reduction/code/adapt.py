from propose import *
started=time.monotonic();labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');Dp=sparse.diags(D)
cs=np.array(json.loads((ROOT/'development/control_proposal.json').read_text())['controls']);cs=np.rint(cs*1e6)/1e6
# A new basis tailored to the selected intervention, not claimed independent of discovery.
sn=[psi[:,None]]
for perturb in [0,-.02,.02]:
 pp=psi.astype(complex)
 for u,v in cs:
  A=K+(u+perturb)*Dp+(v-perturb)*W
  arr=expm_multiply(-.5j*A,pp,start=0,stop=1,num=17,traceA=-.5j*A.diagonal().sum())
  sn.extend([arr.real.T,arr.imag.T]);pp=arr[-1]
X=np.concatenate(sn,axis=1);P,s,_=linalg.svd(X,full_matrices=False,check_finite=False)
np.savez_compressed(ROOT/'development/adapted_basis.npz',J=P,s=s,controls=cs)
print('spectrum',[(d,float(s[d])) for d in [4,8,12,16,24,32,48,64]],flush=True)
for d in [8,12,16,24,32,48,64]:
 J=P[:,:d];hs=np.array([J.T@K@J,J.T@(D[:,None]*J),J.T@W@J]);Es=[K@J-J@hs[0],D[:,None]*J-J@hs[1],W@J-J@hs[2]];zs=[]
 z=J.T@psi;err=np.linalg.norm(psi-J@z)
 for u,v in cs:
  A=hs[0]+u*hs[1]+v*hs[2];ev,V=linalg.eigh(A);grid=np.linspace(0,.5,81)
  arr=(V@(np.exp(-1j*ev[:,None]*grid)*(V.T@z)[:,None])).T;z=arr[-1]
  E=Es[0]+u*Es[1]+v*Es[2];R=E.T@E
  vals=np.maximum(np.einsum('ti,ij,tj->t',arr.conj(),R,arr).real,0)
  err+=np.trapezoid(np.sqrt(vals),grid)
 print('d',d,'D',float(np.real(z.conj()@hs[1]@z)),'err_estimate',err,flush=True)
print('total',time.monotonic()-started)
