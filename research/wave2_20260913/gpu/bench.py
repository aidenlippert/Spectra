import json,time,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
import cupy as cp
import cupyx.scipy.sparse as csp
from cupyx.scipy.sparse.linalg import eigsh as geigsh
sys.path.insert(0,'/tmp/w/src'); from marginal_determinant_tree import DeterminantOracle
R=Path('/home/ubuntu/spectra-wave2-20260913/gpu'); out=R/'results'; out.mkdir(exist_ok=True)
class Unlimited(dict):
 def __len__(self): return 0
for tag in ('h8','h10'):
 c=json.loads((R/'src'/f'{tag}.json').read_text()); w=json.loads((R/'src'/f'{tag}_ref.json').read_text())['independent_upper']; basis=sorted(w['states']); o=DeterminantOracle(c); o.cache=Unlimited(); ix={s:i for i,s in enumerate(basis)}; t=time.perf_counter(); rows=[];cols=[];dat=[]
 for j,s in enumerate(basis):
  for q,v in o.action(s).items():
   if q in ix: rows.append(ix[q]);cols.append(j);dat.append(float(v))
 A=csr_matrix((dat,(rows,cols)),shape=(len(basis),len(basis))); A=(A+A.T)*.5; asm=time.perf_counter()-t
 def cpu(): return eigsh(A,k=1,which='SA',tol=1e-10,maxiter=4000,v0=np.ones(len(basis)))
 D=csp.csr_matrix(A); cp.cuda.Stream.null.synchronize()
 def gpu():
  v,x=geigsh(D,k=1,which='SA',tol=1e-10,maxiter=4000,v0=cp.ones(len(basis))); cp.cuda.Stream.null.synchronize(); return float(v[0].get()),cp.asnumpy(x[:,0])
 def timed(fn):
  z=[]
  for k in range(4):
   t=time.perf_counter(); r=fn();
   if k:z.append(time.perf_counter()-t)
  return r,z
 cr,ct=timed(cpu); gr,gt=timed(gpu); rec={'tag':tag,'n':len(basis),'nnz':A.nnz,'assembly_seconds':asm,'cpu':{'times':ct,'median':float(np.median(ct)),'energy':float(cr[0])},'gpu':{'times':gt,'median':float(np.median(gt)),'energy':float(gr[0])},'energy_absdiff':abs(float(cr[0])-float(gr[0])),'cpu_residual':float(np.linalg.norm(A@cr[1]-cr[0]*cr[1])),'gpu_residual':float(np.linalg.norm(A@gr[1]-gr[0]*gr[1])),'matrix_sha256':hashlib.sha256(np.asarray(A.data).tobytes()).hexdigest()}; (out/f'{tag}.json').write_text(json.dumps(rec,indent=2)); print(rec,flush=True)
