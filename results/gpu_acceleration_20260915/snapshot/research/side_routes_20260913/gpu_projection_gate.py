"""Isolated actual-block CPU/GPU PSD projection benchmark; no solver changes."""
import argparse,hashlib,json,socket,statistics,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import cupy as cp
import cupyx


def run(source,out):
    arrays=np.load(source,allow_pickle=False);cases=[]
    for size in (225,725):
        keys=[k for k in arrays.files if k.startswith('gram_') and arrays[k].shape==(size,size)]
        if len(keys)!=8:raise ValueError('Expected eight actual molecular blocks')
        A=np.stack([(arrays[k]+arrays[k].T)/2 for k in keys])
        reference=[]
        for a in A:
            w,v=eigh(a,driver='evd',check_finite=False);reference.append((v*np.maximum(w,0))@v.T)
        reference=np.stack(reference);methods={}
        for name in ('cpu_evd','gpu_loop','gpu_batch'):
            times=[];result=None;error=None
            try:
                for iteration in range(5):
                    cp.cuda.Stream.null.synchronize();start=time.perf_counter()
                    if name=='cpu_evd':
                        result=[]
                        for a in A:
                            w,v=eigh(a,driver='evd',check_finite=False);result.append((v*np.maximum(w,0))@v.T)
                        result=np.stack(result)
                    else:
                        device=cp.asarray(A)
                        with cupyx.errstate(linalg='raise'):
                            if name=='gpu_loop':
                                projected=[]
                                for a in device:
                                    w,v=cp.linalg.eigh(a);projected.append((v*cp.maximum(w,0))@v.T)
                                result=cp.asnumpy(cp.stack(projected))
                            else:
                                w,v=cp.linalg.eigh(device);projected=(v*cp.maximum(w,0)[:,None,:])@v.swapaxes(-1,-2)
                                result=cp.asnumpy(projected)
                        cp.cuda.Stream.null.synchronize()
                    elapsed=time.perf_counter()-start
                    if iteration>=2:times.append(elapsed)
                difference=float(np.max(np.abs(result-reference)))
                if difference>1e-10:raise AssertionError('PSD projection disagreement')
                methods[name]={'median_seconds':statistics.median(times),'trials_seconds':times,'max_abs_projection_difference':difference,'accepted_kernel_comparison':True}
            except Exception as exc:
                methods[name]={'accepted_kernel_comparison':False,'error':str(exc)}
            print(size,name,methods[name],flush=True)
        cases.append({'size':size,'batch':8,'source_keys':keys,'methods':methods})
        receipt={'host':socket.gethostname(),'device':cp.cuda.runtime.getDeviceProperties(0)['name'].decode(),
                 'cupy':cp.__version__,'numpy':np.__version__,'source_sha256':hashlib.sha256(Path(source).read_bytes()).hexdigest(),
                 'warmups':2,'timed_trials':3,'cases':cases,'gpu_pool_bytes':cp.get_default_memory_pool().total_bytes(),
                 'scope':'Float64 PSD projection only, including reconstruction and GPU input/output transfers. No solver speedup or exact-certificate claim.'}
        Path(out).write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(a.source,a.out)
