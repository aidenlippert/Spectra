"""Bounded grouped A100 screen, with CPU checks on current physical tensors."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from scipy.linalg import eigvalsh
import cupy as cp
import cupyx


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--data',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--batch',type=int,default=256);a=ap.parse_args()
 start=time.perf_counter();data=np.load(a.data,allow_pickle=False);indices=sorted(int(k[1:]) for k in data.files if k.startswith('A'))
 X=data['points_dense'];Y=data['points_diag'];N=len(X)
 if len(indices)!=94 or not 1<=N<=2048 or not 1<=a.batch<=256:raise ValueError('Bounded complete physical screen required')
 grouped={};scalar=[]
 for i in indices:
  A,D,T=data[f'A{i}'],data[f'D{i}'],data[f'T{i}'];scalar.append((A,D,T));grouped.setdefault(len(A),[]).append((A,D,T))
 groups=[];sync=cp.cuda.Stream.null.synchronize
 for n,items in sorted(grouped.items()):
  base=np.stack([v[0] for v in items]);D=np.stack([v[1] for v in items]).transpose(1,0,2,3).reshape(8,-1);T=np.stack([v[2] for v in items]).transpose(1,0,2).reshape(Y.shape[1],-1)
  groups.append((n,cp.asarray(base),cp.asarray(D),cp.asarray(T)))
 dX,dY=cp.asarray(X),cp.asarray(Y);sync();setup=time.perf_counter()-start
 def evaluate(x,y):
  lowest=cp.full(len(x),cp.inf)
  for n,A,D,T in groups:
   M=(x@D).reshape(len(x),len(A),n,n)+A
   diag=(y@T).reshape(len(x),len(A),n);idx=cp.arange(n);M[...,idx,idx]+=diag
   vals=M[...,0,0] if n==1 else cp.linalg.eigvalsh(M.reshape(-1,n,n))[:,0].reshape(len(x),len(A))
   lowest=cp.minimum(lowest,vals.min(axis=1))
  return lowest
 warm=time.perf_counter()
 with cupyx.errstate(linalg='raise'):
  evaluate(dX[:min(N,a.batch)],dY[:min(N,a.batch)]);sync()
 warm=time.perf_counter()-warm;run=time.perf_counter()
 with cupyx.errstate(linalg='raise'):
  got=cp.asnumpy(cp.concatenate([evaluate(dX[j:j+a.batch],dY[j:j+a.batch]) for j in range(0,N,a.batch)]));sync()
 gpu_seconds=time.perf_counter()-run
 costs=(data['base_penalties'][None]+X[:,6:8])@data['theta'];densities=(got-costs)/5;order=np.argsort(-densities)
 sample=sorted(set([0,N-1,*order[:12].tolist(),*np.random.default_rng(7).choice(N,min(12,N),replace=False).tolist()]))
 cpu_start=time.perf_counter();reference=[]
 for j in sample:
  reference.append(min(float(eigvalsh(A+np.einsum('i,ijk->jk',X[j],D)+np.diag(Y[j]@T),subset_by_index=[0,0],check_finite=False)[0]) for A,D,T in scalar))
 error=float(np.max(abs(got[sample]-reference)))
 if not np.allclose(got[sample],reference,atol=1e-9,rtol=1e-10):raise ValueError('GPU/CPU current-model mismatch')
 result={'numerical_only':True,'data_sha256':hashlib.sha256(a.data.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'candidates':N,'sectors':94,'group_sizes':len(groups),'setup_transfer_seconds':setup,'warmup_seconds':warm,'gpu_seconds':gpu_seconds,'cpu_check_count':len(sample),'cpu_check_seconds':time.perf_counter()-cpu_start,'cpu_max_abs_error':error,'baseline_density':float(densities[0]),'best_density':float(densities[order[0]]),'best_index':int(order[0]),'best_minimum':float(got[order[0]]),'top_indices':order[:20].tolist(),'density_values':densities.tolist(),'minimum_values':got.tolist(),'resident_bytes':int(sum(A.nbytes+D.nbytes+T.nbytes for _,A,D,T in groups)+dX.nbytes+dY.nbytes),'max_matrix_batch_bytes':max(a.batch*len(A)*n*n*8 for n,A,D,T in groups),'device_name':str(cp.cuda.runtime.getDeviceProperties(0)['name']),'scope':'All94 sectors per candidate; LP-directed extra-overlap screening. CPU checks cover winners, baseline, boundary and random candidates. No exact positivity or ground-energy acceptance.'}
 a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('density_values','minimum_values')}),flush=True)


if __name__=='__main__':main()
