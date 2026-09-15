import json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh
import cupy as cp
import cupyx.scipy.sparse as csp
from cupyx.scipy.sparse.linalg import eigsh as geigsh
from research.all_angles_20260913.selected_refinement.refine import CompiledHamiltonian,restricted_matrix
R=Path('/home/ubuntu/spectra-wave2-20260913/gpu'); f=R/'src/h10.json'; fixture=json.loads(f.read_text()); w=json.loads((R/'src/h10_ref.json').read_text())['independent_upper']; basis=sorted(w['states']); t=time.perf_counter(); A=restricted_matrix(CompiledHamiltonian(fixture),basis); assembly=time.perf_counter()-t
np.savez(R/'results/h10_csr.npz',data=A.data,indices=A.indices,indptr=A.indptr,shape=A.shape); sh=hashlib.sha256((R/'results/h10_csr.npz').read_bytes()).hexdigest()
def cpu(): return eigsh(A,k=1,which='SA',tol=1e-10,maxiter=4000,v0=np.ones(len(basis)))
D=csp.csr_matrix(A); cp.cuda.Stream.null.synchronize()
transfer=[]
for _ in range(3):
 t0=time.perf_counter(); _d=csp.csr_matrix(A); cp.cuda.Stream.null.synchronize(); transfer.append(time.perf_counter()-t0)
def gpu():
 v,x=geigsh(D,k=1,which='SA',tol=1e-10,maxiter=4000,v0=cp.ones(len(basis))); cp.cuda.Stream.null.synchronize(); return float(v[0].get()),cp.asnumpy(x[:,0])
def tm(fn):
 z=[]
 for k in range(4):
  t=time.perf_counter(); r=fn()
  if k:z.append(time.perf_counter()-t)
 return r,z
cr,ct=tm(cpu); gr,gt=tm(gpu); amps=[int(round(float(x)*10**12)) for x in gr[1]]; witness={'states':[s for s,a in zip(basis,amps) if a],'amplitudes':[a for a in amps if a]}; from research.certificate_scaling.streaming_reference_upper import upper; exact,stats=upper(fixture,witness); (R/'results/h10_gpu_witness.json').write_text(json.dumps({'independent_upper':witness,'upper':str(exact)},separators=(',',':'))+'\n'); rec={'n':len(basis),'nnz':A.nnz,'assembly_seconds':assembly,'transfer_times':transfer,'cpu_times':ct,'gpu_times':gt,'cpu_median':float(np.median(ct)),'gpu_median':float(np.median(gt)),'cpu_energy':float(cr[0]),'gpu_energy':float(gr[0]),'energy_absdiff':abs(float(cr[0])-float(gr[0])),'cpu_residual':float(np.linalg.norm(A@cr[1]-cr[0]*cr[1])),'gpu_residual':float(np.linalg.norm(A@gr[1]-gr[0]*gr[1])),'csr_sha256':sh,'witness_sha256':hashlib.sha256((R/'results/h10_gpu_witness.json').read_bytes()).hexdigest(),'upper':str(exact),'upper_stats':stats}; (R/'results/h10_receipt.json').write_text(json.dumps(rec,indent=2)+'\n'); print(rec)
