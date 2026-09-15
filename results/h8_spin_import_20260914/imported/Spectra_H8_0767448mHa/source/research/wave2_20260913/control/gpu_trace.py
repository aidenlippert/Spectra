import json, time, hashlib
from pathlib import Path
import numpy as np
import scipy.sparse as sp

KERNEL = r'''extern "C" __global__ void emit(const unsigned long long* states,const unsigned long long* req,const unsigned long long* occ,const unsigned long long* flip,const unsigned long long* parity,const long long* coef,int n,int nt,const unsigned long long* basis,int* rows,long long* vals){int k=blockDim.x*blockIdx.x+threadIdx.x; long long total=(long long)n*nt; if(k>=total)return; int j=k/nt,t=k%nt; unsigned long long s=states[j]; rows[k]=-1; if((s&req[t])!=occ[t])return; unsigned long long x=s^flip[t]; int lo=0,hi=n; while(lo<hi){int m=(lo+hi)/2; if(basis[m]<x)lo=m+1;else hi=m;} if(lo<n&&basis[lo]==x){int p=__popcll(s&parity[t]); rows[k]=lo; vals[k]=(p&1)?-coef[t]:coef[t];}}'''

def run(fixture_path, basis_path, cpu_path, out_path):
    import cupy as cp, cupyx.scipy.sparse as csp
    from research.all_angles_20260913.selected_refinement.refine import CompiledHamiltonian
    fixture_raw=Path(fixture_path).read_bytes(); fixture=json.loads(fixture_raw)
    basis=sorted(json.loads(Path(basis_path).read_text())['independent_upper']['states']); n=len(basis)
    comp=CompiledHamiltonian(fixture); D=comp.denominator; terms=comp.terms; nt=len(terms)
    if sum(abs(t[4]) for t in terms)>=2**53: raise RuntimeError('integer coefficient bound refused')
    a=np.asarray(basis,dtype=np.uint64); req=np.array([t[0] for t in terms],dtype=np.uint64); occ=np.array([t[1] for t in terms],dtype=np.uint64); flip=np.array([t[2] for t in terms],dtype=np.uint64); parity=np.array([t[3] for t in terms],dtype=np.uint64); coef=np.array([t[4] for t in terms],dtype=np.int64)
    mod=cp.RawModule(code=KERNEL,options=('-std=c++11',)); ker=mod.get_function('emit')
    def assemble():
        total=n*nt; drows=cp.full(total,-1,dtype=cp.int32); dvals=cp.zeros(total,dtype=cp.int64); ds=cp.asarray(a); db=ds
        ker(((total+255)//256,),(256,),(ds,cp.asarray(req),cp.asarray(occ),cp.asarray(flip),cp.asarray(parity),cp.asarray(coef),n,nt,db,drows,dvals))
        print("KERNEL_DIAG",int(cp.sum(dvals[:nt][drows[:nt]==0]).get()),"CPU_EXPECT",int(round(comp.diagonal(basis[0])*D)),flush=True)
        mask=drows>=0; rows=cp.asnumpy(drows[mask]); vals=cp.asnumpy(dvals[mask]); cols=np.repeat(np.arange(n,dtype=np.int32),nt)[cp.asnumpy(mask)];
        # CuPy sparse accepts floating payloads; the validated <2**53 bound makes
        # these integer numerators exactly representable in float64.
        M = csp.coo_matrix((cp.asarray(vals,dtype=cp.float64),(cp.asarray(rows),cp.asarray(cols))),shape=(n,n)).tocsr()
        print("CSR_BEFORE",M.diagonal()[:3].get().tolist(),"NNZ",M.nnz,flush=True)
        M.sum_duplicates(); M.eliminate_zeros(); M.data /= D
        print("CSR_AFTER",M.diagonal()[:3].get().tolist(),"NNZ",M.nnz,flush=True)
        return M
    t=time.perf_counter(); A=assemble(); cp.cuda.Stream.null.synchronize(); cold=time.perf_counter()-t
    warm=[]
    for _ in range(0): t=time.perf_counter(); B=assemble(); cp.cuda.Stream.null.synchronize(); warm.append(time.perf_counter()-t)
    ref=np.load(cpu_path); cpu=sp.csr_matrix((ref['data'],ref['indices'],ref['indptr']),shape=tuple(ref['shape']))
    host = A.get()
    host.sum_duplicates(); host.eliminate_zeros(); host.sort_indices()
    cpu.sum_duplicates(); cpu.eliminate_zeros(); cpu.sort_indices()
    difference = host-cpu
    difference.eliminate_zeros()
    print("DIAGNOSTIC", json.dumps({"shape_match":host.shape==cpu.shape,"indptr_match":np.array_equal(host.indptr,cpu.indptr),"indices_match":np.array_equal(host.indices,cpu.indices),"operator_max_abs_error":float(np.max(np.abs(difference.data),initial=0)),"difference_nnz":difference.nnz,"gpu_diagonal":host.diagonal()[:5].tolist(),"cpu_diagonal":cpu.diagonal()[:5].tolist(),"gpu_row0_indices":host.getrow(0).indices[:8].tolist(),"cpu_row0_indices":cpu.getrow(0).indices[:8].tolist()}),flush=True)
    sp.save_npz(str(out_path)+".gpu.npz",host)
    maxdiff=float(np.max(np.abs(host.data-cpu.data))) if host.nnz == cpu.nnz else float("inf")
    eq=(A.shape==cpu.shape and np.array_equal(A.indptr.get(),cpu.indptr) and np.array_equal(A.indices.get(),cpu.indices) and maxdiff <= np.finfo(np.float64).eps)
    rec={'n':n,'terms':nt,'denominator':D,'sum_abs_integer_coeff':int(sum(abs(t[4]) for t in terms)),'matrix_nnz':int(A.nnz),'shape':A.shape,'cold_jit_assemble_seconds':cold,'warm_assemble_seconds':warm,'gpu_memory_bytes_final_snapshot':int(cp.cuda.Device().mem_info[1]-cp.cuda.Device().mem_info[0]),'host_transfer_bytes':int(A.data.nbytes+A.indices.nbytes+A.indptr.nbytes),'cpu_csr_match':bool(eq),'max_abs_data_error':maxdiff,'gpu_csr_sha256':hashlib.sha256(A.data.get().tobytes()+A.indices.get().tobytes()+A.indptr.get().tobytes()).hexdigest(),'cpu_csr_sha256':hashlib.sha256(cpu.data.tobytes()+cpu.indices.tobytes()+cpu.indptr.tobytes()).hexdigest()}
    Path(out_path).write_text(json.dumps(rec,indent=2)+'\n'); return rec

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(); p.add_argument('--fixture',required=True);p.add_argument('--basis',required=True);p.add_argument('--cpu',required=True);p.add_argument('--out',required=True);a=p.parse_args(); print(run(a.fixture,a.basis,a.cpu,a.out))
