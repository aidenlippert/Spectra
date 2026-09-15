"""Isolated dsyev/dsyevd gate using SCS's actual bundled OpenBLAS."""
from pathlib import Path
import argparse,ctypes as C,hashlib,json,time,statistics,socket
import numpy as np
import scs,scs._scs_direct


def run(raw_path,out):
    library=next((Path(scs.__file__).parent.parent/'scs.libs').glob('libopenblas*.so'))
    lib=C.CDLL(str(library));lib.openblas_get_config.restype=C.c_char_p
    config=lib.openblas_get_config().decode()
    if 'USE64BITINT' in config:raise ValueError('This harness requires LP64 LAPACK')
    lib.openblas_set_num_threads.argtypes=[C.c_int];lib.openblas_set_num_threads(1)
    lib.openblas_get_num_threads.restype=C.c_int
    if lib.openblas_get_num_threads()!=1:raise AssertionError('Thread count not applied')
    ptr=lambda a:a.ctypes.data_as(C.POINTER(C.c_double))
    iptr=lambda a:a.ctypes.data_as(C.POINTER(C.c_int))
    def workspace(n,driver):
        a=np.eye(n,dtype=np.float64,order='F');w=np.empty(n);work=np.empty(1);iw=np.empty(1,dtype=np.int32)
        nn=C.c_int(n);lda=C.c_int(n);lw=C.c_int(-1);liw=C.c_int(-1);info=C.c_int()
        args=[C.c_char_p(b'V'),C.c_char_p(b'L'),C.byref(nn),ptr(a),C.byref(lda),ptr(w),ptr(work),C.byref(lw)]
        if driver=='dsyevd_':args += [iptr(iw),C.byref(liw)]
        args += [C.byref(info)];getattr(lib,driver)(*args)
        if info.value:raise RuntimeError(('workspace',driver,info.value))
        return np.empty(int(work[0])),np.empty(int(iw[0]) if driver=='dsyevd_' else 1,dtype=np.int32)
    def call(a,driver,work,iw):
        n=C.c_int(len(a));lda=C.c_int(len(a));w=np.empty(len(a));lw=C.c_int(len(work));liw=C.c_int(len(iw));info=C.c_int()
        args=[C.c_char_p(b'V'),C.c_char_p(b'L'),C.byref(n),ptr(a),C.byref(lda),ptr(w),ptr(work),C.byref(lw)]
        if driver=='dsyevd_':args += [iptr(iw),C.byref(liw)]
        args += [C.byref(info)];getattr(lib,driver)(*args)
        if info.value:raise RuntimeError(('solve',driver,info.value))
        return w,a
    for driver in ('dsyev_','dsyevd_'):
        work,iw=workspace(2,driver);w,u=call(np.array([[2.,1.],[1.,2.]],order='F'),driver,work,iw)
        assert np.max(abs(w-np.array([1.,3.])))<1e-12
    arrays=np.load(raw_path,allow_pickle=False);cases=[]
    for size in (725,225):
        names=[k for k in arrays.files if k.startswith('gram_') and arrays[k].shape==(size,size)]
        if len(names)!=8:raise ValueError(('Expected eight actual blocks',size,names))
        originals=[(arrays[k]+arrays[k].T)/2 for k in names];results={};projections={};errors={}
        for driver in ('dsyev_','dsyevd_'):
            work,iw=workspace(size,driver);times=[]
            for trial in range(5):
                inputs=[np.array(a,order='F',copy=True) for a in originals];start=time.perf_counter()
                eigenpairs=[call(a,driver,work,iw) for a in inputs]
                elapsed=time.perf_counter()-start
                if trial>=2:times.append(elapsed)
            err=max(float(np.max(np.abs((u*w)@u.T-a))) for a,(w,u) in zip(originals,eigenpairs))
            if err>1e-10:raise AssertionError(('Eigen residual',driver,err))
            errors[driver]=err;projections[driver]=[(u*np.maximum(w,0))@u.T for w,u in eigenpairs]
            results[driver]={'seconds':times,'median_seconds':statistics.median(times)}
        agreement=max(float(np.max(np.abs(a-b))) for a,b in zip(projections['dsyev_'],projections['dsyevd_']))
        if agreement>1e-10:raise AssertionError(('PSD projection disagreement',agreement))
        cases.append({'batch':8,'size':size,'source_keys':names,'drivers':results,'eigen_reconstruction_error':errors,'PSD_projection_disagreement':agreement})
        print(json.dumps(cases[-1]),flush=True)
    receipt={'host':socket.gethostname(),'scs':scs.__version__,'library':str(library),'library_sha256':hashlib.sha256(library.read_bytes()).hexdigest(),'openblas_config':config,'threads':lib.openblas_get_num_threads(),'source':str(raw_path),'source_sha256':hashlib.sha256(Path(raw_path).read_bytes()).hexdigest(),'warmups':2,'timed_trials':3,'cases':cases,'scope':'Isolated actual-block eigenkernel timing; copies, workspace allocation, reconstruction and complete solver overhead excluded from timed kernel.'}
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.raw,a.out)
