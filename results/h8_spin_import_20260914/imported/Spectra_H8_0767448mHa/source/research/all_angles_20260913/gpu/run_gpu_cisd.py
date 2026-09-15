import sys,json,time,hashlib,socket
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0,str(Path(__file__).parent/'src'))
from marginal_determinant_tree import DeterminantOracle
import cupy as cp

ROOT=Path(__file__).parent
def one(tag):
 f=ROOT/'src'/f'{tag}.json'; cert=json.loads(f.read_text()); oracle=DeterminantOracle(cert)
 states=[s for s in range(1<<oracle.modes) if oracle.valid_state(s)]
 def de(s): return sum(float(c) for support,c in oracle.diagonal.items() if s&support==support)
 hf=min(states,key=de); occ=[i for i in range(oracle.modes) if hf>>i&1]; vir=[i for i in range(oracle.modes) if not hf>>i&1]; basis={hf}
 for i in occ:
  for a in vir:basis.add(hf^(1<<i)^(1<<a))
 for x,i in enumerate(occ):
  for j in occ[x+1:]:
   for y,a in enumerate(vir):
    for b in vir[y+1:]:basis.add(hf^(1<<i)^(1<<j)^(1<<a)^(1<<b))
 basis=sorted(basis); idx={s:i for i,s in enumerate(basis)}; n=len(basis); A=np.zeros((n,n),dtype=np.float64); t=time.perf_counter()
 for j,s in enumerate(basis):
  for target,v in oracle.action(s).items():
   if target in idx:A[idx[target],j]=float(v)
 A=(A+A.T)/2; assemble=time.perf_counter()-t; np.save(ROOT/'results'/f'{tag}_H.npy',A)
 sha=hashlib.sha256((ROOT/'results'/f'{tag}_H.npy').read_bytes()).hexdigest(); out={"tag":tag,"n":n,"assemble_seconds":assemble,"matrix_sha256":sha,"basis_sha256":hashlib.sha256(np.asarray(basis,dtype=np.uint64).tobytes()).hexdigest()}
 # CPU baseline and GPU eigh; 3 warmups + 3 reps, transfer included
 for name in ['cpu','gpu']:
  ts=[]; vals=None
  for k in range(6):
   t=time.perf_counter()
   if name=='cpu': vals=eigh(A,driver='evr',check_finite=False,subset_by_index=[0,0],eigvals_only=False)
   else:
    d=cp.asarray(A); cp.cuda.Stream.null.synchronize(); w,v=cp.linalg.eigh(d); x=float(w[0].get()); vec=cp.asnumpy(v[:,0]); cp.cuda.Stream.null.synchronize(); vals=(x,vec)
   dt=time.perf_counter()-t
   if k>=3:ts.append(dt)
  out[name]={"times":ts,"median_seconds":float(np.median(ts)),"energy":float(vals[0])}
 # exact rational streaming replay of rationalized GPU vector
 scale=10**12; amps=[int(round(float(x)*scale)) for x in (vec if 'vec' in locals() else vals[1])]; wit={"states":[s for s,a in zip(basis,amps) if a],"amplitudes":[a for a in amps if a]}
 from streaming_reference_upper import upper
 val,stats=upper(cert,wit); out.update({"upper_rational":str(val),"upper_stats":stats,"witness_support":len(wit['states']),"gpu_cpu_energy_absdiff":abs(out['gpu']['energy']-out['cpu']['energy'])})
 np.savez(ROOT/'results'/f'{tag}_gpu_witness.npz',states=np.asarray(wit['states'],dtype=np.uint64),amplitudes=np.asarray(wit['amplitudes'],dtype=np.int64),basis=np.asarray(basis,dtype=np.uint64))
 (ROOT/'results'/f'{tag}_receipt.json').write_text(json.dumps(out,indent=2)+'\n'); return out
if __name__=='__main__':
 rows=[]
 for tag in ('h8','h10'):
  try: rows.append(one(tag))
  except Exception as e: rows.append({'tag':tag,'error':repr(e)})
 (ROOT/'results'/'summary.json').write_text(json.dumps(rows,indent=2)+'\n'); print(json.dumps(rows,indent=2))
