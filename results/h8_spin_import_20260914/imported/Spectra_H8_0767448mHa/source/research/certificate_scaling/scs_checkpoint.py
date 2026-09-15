"""Bounded SCS iterate restarts with canonical-data identity checks.

This saves even optimal_inaccurate finite iterates. It does not preserve
SCS's internal scaling/acceleration workspace or certify solver convergence.
"""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
import cvxpy as cp
import scs
from cvxpy.reductions.solvers.conic_solvers.scs_conif import dims_to_solver_dict


def array_hash(value):
    a=np.ascontiguousarray(value)
    h=hashlib.sha256();h.update(str(a.dtype).encode());h.update(json.dumps(a.shape).encode());h.update(a.tobytes())
    return h.hexdigest()


def identity(data,options,contract):
    matrices={}
    for key in ('A','P'):
        if key in data and data[key] is not None:
            a=data[key].tocsc(copy=True);a.sum_duplicates();a.sort_indices()
            matrices[key]={'shape':list(a.shape),'data':array_hash(a.data),'indices':array_hash(a.indices.astype(np.int64)),'indptr':array_hash(a.indptr.astype(np.int64))}
    cones=dims_to_solver_dict(data['dims'])
    record={'format':'scs_canonical_restart_v1','cvxpy':cp.__version__,'scs':scs.__version__,'numpy':np.__version__,
            'matrices':matrices,'b':array_hash(data['b']),'c':array_hash(data['c']),'cones':cones,
            'options':{k:v for k,v in options.items() if k not in ('time_limit_secs','max_iters')},
            'contract':contract}
    # An isolated numerical-backend experiment must not silently resume a
    # stock checkpoint just because the Python package version is unchanged.
    overrides={}
    for variable in ('LD_PRELOAD','SPECTRA_EVD_LIBRARY'):
        if os.environ.get(variable):
            entries=[]
            for name in os.environ[variable].replace(':',' ').split():
                path=Path(name)
                if not path.is_absolute() or not path.is_file():
                    raise ValueError('Checkpoint backend overrides require explicit library paths')
                entries.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            overrides[variable]=entries
    if overrides:record['declared_backend_overrides']=overrides
    # JSON roundtrip normalizes tuples and NumPy-free scalar representations.
    return json.loads(json.dumps(record,sort_keys=True))


def state_paths(prefix):
    prefix=Path(prefix);return Path(str(prefix)+'.json'),Path(str(prefix)+'.npz')


def load(prefix,expected,data):
    meta_path,arrays_path=state_paths(prefix);meta=json.loads(meta_path.read_text())
    if meta['identity']!=expected:raise ValueError('SCS checkpoint canonical identity/version/options mismatch')
    if hashlib.sha256(arrays_path.read_bytes()).hexdigest()!=meta['arrays_sha256']:raise ValueError('SCS checkpoint array hash mismatch')
    with np.load(arrays_path,allow_pickle=False) as stored:
        if set(stored.files)!={'x','y','s'}:raise ValueError('Unexpected checkpoint array keys')
        result={key:stored[key].copy() for key in ('x','y','s')}
    n=len(data['c']);m=len(data['b'])
    for key,size in [('x',n),('y',m),('s',m)]:
        if result[key].shape!=(size,) or result[key].dtype!=np.float64 or not np.all(np.isfinite(result[key])):
            raise ValueError('Invalid finite canonical SCS iterate')
    if meta['status_val'] not in (1,2):raise ValueError('Checkpoint is not a primal/dual solution iterate')
    return result,meta


def solve(problem,options,output_prefix,resume=None,contract=None):
    start=time.monotonic();meta_path,arrays_path=state_paths(output_prefix)
    if meta_path.exists() or arrays_path.exists():raise FileExistsError('Preserve existing SCS checkpoint; use a new output path')
    data,chain,inverse=problem.get_problem_data(cp.SCS)
    expected=identity(data,options,contract);cache={};parent=None
    if resume is not None:
        state,parent=load(resume,expected,data);cache['SCS']=state
    compiled=time.monotonic()
    result=chain.solver.solve_via_data(data,warm_start=resume is not None,verbose=False,solver_opts=dict(options),solver_cache=cache)
    solved=time.monotonic();problem.unpack_results(result,chain,inverse)
    info=result['info'];status=int(info['status_val'])
    if status not in (1,2) or any(not np.all(np.isfinite(result[k])) for k in ('x','y','s')):
        raise RuntimeError('SCS did not return a finite solution iterate to checkpoint')
    meta_path.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(str(arrays_path)+'.tmp.npz')
    np.savez_compressed(temp,**{k:np.asarray(result[k],dtype=np.float64) for k in ('x','y','s')});os.replace(temp,arrays_path)
    meta={'identity':expected,'status_val':status,'status':str(info['status']),'arrays_sha256':hashlib.sha256(arrays_path.read_bytes()).hexdigest(),
          'resumed_from':str(resume) if resume is not None else None,'iterations_this_stage':int(info['iter']),
          'cumulative_iterations':int(info['iter'])+(parent['cumulative_iterations'] if parent else 0),
          'compile_and_load_seconds':compiled-start,'solve_seconds':solved-compiled,'wall_seconds':time.monotonic()-start,
          'info':{k:(v.item() if isinstance(v,np.generic) else v) for k,v in info.items()},
          'scope':'Canonical x/y/s iterate restart; internal SCS scaling, acceleration and factorization are rebuilt; no certified optimum claim.'}
    temp_meta=Path(str(meta_path)+'.tmp');temp_meta.write_text(json.dumps(meta,indent=2)+'\n');os.replace(temp_meta,meta_path)
    return meta
