from common import *
from scipy.sparse.linalg import lsqr
start=time.monotonic();op=load_operator();T=sparse.load_npz(PREP/'twirl.npz');selected=np.load(PREP/'selected.npy');p=np.load(OUT/'mps_six_dual_coordinates.npy');target=T.T@p
sol=lsqr(T[selected].T,target,atol=1e-13,btol=1e-13,iter_lim=500)
y=sol[0]/op.scale
print('moment_pullback_residual',sol[3], 'iters',sol[2],flush=True)
np.save(OUT/'physical_moment_dual.npy',y)
records=[];Cs={}
for k,(M,V) in enumerate(zip(op.M,op.V)):
 n=len(V);C=np.asarray(M.T@y).reshape(n,n);C=(C+C.T)/2;Cs[f'C_{k}']=C
 ev,E=linalg.eigh(C,check_finite=False)
 s=V.shape[1];base=linalg.orth(V);low=E[:,:min(n,s)]
 residual=low-base@(base.T@low)
 rec={'k':k,'dimension':n,'retained':s,'min_eigenvalue':float(ev[0]),'max_eigenvalue':float(ev[-1]),'near_null':int(np.sum(ev<1e-9)),'low_eigenvalues':ev[[j for j in [0,15,31,63,95,127,159,191,255,n-1] if j<n]].tolist(),'squared_low_space_missing':float(np.sum(residual**2))}
 records.append(rec)
 if n>64:print(rec,flush=True)
 if ev[0]<-1e-8:raise AssertionError('Computed physical moment matrix is not PSD')
np.savez_compressed(OUT/'physical_moment_matrices.npz',**Cs)
record_json(OUT/'physical_moment_probe.json',{'seconds':time.monotonic()-start,'pullback_residual':sol[3],'blocks':records,'guide_only':True})
