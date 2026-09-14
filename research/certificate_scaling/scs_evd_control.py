"""Bounded stock/preload SCS projection controls; JSON plus shim stderr counts."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time
import numpy as np
import cvxpy as cp
import scipy.linalg
results=[]
for n in (4,20):
    rng=np.random.default_rng(n);a=rng.normal(size=(n,n));a=(a+a.T)/2
    values,vectors=np.linalg.eigh(a);expected=float(np.maximum(values,0).sum())
    x=cp.Variable((n,n),PSD=True)
    problem=cp.Problem(cp.Minimize(cp.trace(x)),[x-a>>0])
    start=time.perf_counter();problem.solve(solver='SCS',eps=1e-8,max_iters=10000,time_limit_secs=30)
    error=abs(float(problem.value)-expected)
    feasibility=min(float(np.linalg.eigvalsh(x.value).min()),float(np.linalg.eigvalsh(x.value-a).min()))
    if problem.status!='optimal' or error>1e-5 or feasibility < -1e-6:raise AssertionError((n,problem.status,error,feasibility))
    results.append({'n':n,'status':problem.status,'objective':float(problem.value),'expected':expected,
                    'absolute_error':error,'min_cone_eigenvalue':feasibility,'iterations':problem.solver_stats.num_iters,'wall_seconds':time.perf_counter()-start})
print(json.dumps({'preload':os.environ.get('LD_PRELOAD'),'controls':results},indent=2),flush=True)
