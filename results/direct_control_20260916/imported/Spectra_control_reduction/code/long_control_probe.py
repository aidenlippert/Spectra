from propose import *
labels,H,D,W,psi=setup();K=H+9.255*sparse.eye(H.shape[0],format='csr');Dp=sparse.diags(D); start=time.monotonic()
for u in [-.3,-.2,-.1,0,.1,.2,.3]:
 A=K+u*Dp+.1*W;arr=expm_multiply(-1j*A,psi,start=0,stop=30,num=61,traceA=-1j*A.diagonal().sum());val=(abs(arr)**2)@D
 k=int(np.argmin(val));print(u,'minD',val[k],'time',.5*k,'final',val[-1],flush=True)
print('cost',time.monotonic()-start)
