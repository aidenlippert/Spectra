"""Controlled actual-block SciPy eigenkernel benchmark, not a solver claim."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse, ctypes as C, hashlib, json, socket, statistics, time
from pathlib import Path
import numpy as np
import scipy, scipy.linalg as la


def run(raw,out):
    library=next((Path(scipy.__file__).parent.parent/'scipy.libs').glob('libscipy_openblas*.so'))
    lib=C.CDLL(str(library));lib.scipy_openblas_get_num_threads.restype=C.c_int
    lib.scipy_openblas_get_config.restype=C.c_char_p
    threads=lib.scipy_openblas_get_num_threads()
    if threads!=1:raise AssertionError(('Runtime thread count',threads))
    source=np.load(raw,allow_pickle=False);cases=[]
    for n in (725,225):
        keys=[k for k in source.files if k.startswith('gram_') and source[k].shape==(n,n)]
        if len(keys)!=8:raise ValueError('Expected eight blocks')
        originals=[(source[k]+source[k].T)/2 for k in keys];drivers={};projected={}
        for driver in ('ev','evd'):
            times=[]
            for trial in range(5):
                inputs=[np.array(a,order='F',copy=True) for a in originals]
                start=time.perf_counter()
                pairs=[la.eigh(a,driver=driver,check_finite=False,overwrite_a=True) for a in inputs]
                elapsed=time.perf_counter()-start
                if trial>=2:times.append(elapsed)
            residual=max(float(np.max(abs((v*w)@v.T-a))) for a,(w,v) in zip(originals,pairs))
            if residual>1e-10:raise AssertionError(('Eigen reconstruction',driver,residual))
            projected[driver]=[(v*np.maximum(w,0))@v.T for w,v in pairs]
            drivers[driver]={'eight_block_seconds':times,'batch_median_seconds':statistics.median(times),'max_reconstruction_residual':residual}
        agreement=max(float(np.max(abs(a-b))) for a,b in zip(projected['ev'],projected['evd']))
        if agreement>1e-10:raise AssertionError(('Projection disagreement',agreement))
        cases.append({'batch':8,'size':n,'keys':keys,'drivers':drivers,'PSD_projection_disagreement':agreement})
    receipt={'host':socket.gethostname(),'scipy':scipy.__version__,'library':str(library),
             'library_sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
             'openblas_config':lib.scipy_openblas_get_config().decode(),'verified_threads':threads,
             'source':str(raw),'source_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
             'warmups':2,'timed_trials':3,'cases':cases,
             'scope':'Eight-block eigenkernel batches; excludes input copying and reconstruction, includes SciPy call/workspace overhead. No full solver claim.'}
    out.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.raw,a.out)
