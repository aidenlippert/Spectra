import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from propose import *
from scipy.optimize import minimize
started=time.monotonic();labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');Dp=sparse.diags(D)
sn=[psi[:,None]];rng=np.random.default_rng(9135);paths=[]
for k in range(16):
 if k<4: cs=np.tile([[-.1 if k//2 else .1, -.5 if k%2 else .5]],(4,1))
 else:cs=np.column_stack([rng.uniform(-.1,.1,4),rng.uniform(-.5,.5,4)])
 paths.append(cs.tolist());pp=psi.astype(complex)
 for u,v in cs:
  A=K+u*Dp+v*W;arr=expm_multiply(-.5j*A,pp,start=0,stop=1,num=9,traceA=-.5j*A.diagonal().sum())
  sn.extend([arr.real.T,arr.imag.T]);pp=arr[-1]
X=np.concatenate(sn,axis=1);P,s,_=linalg.svd(X,full_matrices=False,check_finite=False)
np.savez_compressed(ROOT/'development'/'short_basis_proposal.npz',J=P,s=s)
print('Svals',[(r,float(s[r])) for r in [8,16,24,32,48,64,96,128,160]],'seconds',time.monotonic()-started,flush=True)
cs=np.array([[.07,.45],[-.09,.48],[.09,.5],[-.03,.42]])
for d in [8,16,24,32,48,64,96]:
 J=P[:,:d];h=[J.T@K@J,J.T@(D[:,None]*J),J.T@W@J];Es=[K@J-J@h[0],D[:,None]*J-J@h[1],W@J-J@h[2]]
 def evalcs(cs,check=False):
  z=J.T@psi;pp=psi.astype(complex);eta=np.linalg.norm(psi-J@z)
  for u,v in cs:
   hu=h[0]+u*h[1]+v*h[2];ev,V=linalg.eigh(hu);grid=np.linspace(0,.5,51)
   traj=(V@(np.exp(-1j*ev[:,None]*grid)*(V.T@z)[:,None])).T;z=traj[-1]
   E=Es[0]+u*Es[1]+v*Es[2];Dg=E.T@E
   res=np.sqrt(np.maximum(np.einsum('ti,ij,tj->t',traj.conj(),Dg,traj).real,0));eta+=np.trapezoid(res,grid)
   if check:
    AA=K+u*Dp+v*W;pp=expm_multiply(-.5j*AA,pp,traceA=-.5j*AA.diagonal().sum())
  return float((z.conj()@h[1]@z).real),eta,np.linalg.norm(pp-J@z) if check else 0
 v,e,err=evalcs(cs,True);print('test d',d,'D',v,'errboundestimate',e,'actual',err,flush=True)
 if d==32:
  def obj(x):
   z=J.T@psi
   for u,v in x.reshape(4,2):
    ev,V=linalg.eigh(h[0]+u*h[1]+v*h[2]);z=V@(np.exp(-.5j*ev)*(V.T@z))
   return float((z.conj()@h[1]@z).real)
  sol=minimize(obj,cs.ravel(),method='L-BFGS-B',bounds=[(-.5,.5),(-.5,.5)]*4,options={'maxiter':100,'ftol':1e-13})
  cs=sol.x.reshape(4,2);print('optimized',cs.tolist(),evalcs(cs,True),sol.message,flush=True)
  (ROOT/'development'/'control_proposal.json').write_text(json.dumps({'controls':cs.tolist(),'training':paths,'seconds':time.monotonic()-started},indent=2))
print('total',time.monotonic()-started)
