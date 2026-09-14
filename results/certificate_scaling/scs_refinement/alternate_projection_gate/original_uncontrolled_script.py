import numpy as np,time,statistics,json,scipy,scipy.linalg as la,os
os.environ["OPENBLAS_NUM_THREADS"]="1";z=np.load("/home/ubuntu/spectra-precision-v2/results/h10_scs_quotient_conditioned/raw_round_0.npz");out=[]
for n in (725,225):
 xs=[z[f"gram_{i}"] for i in range(50) if z[f"gram_{i}"].shape==(n,n)][:8]
 for d in ("ev","evd"):
  ts=[];rr=[]
  for a in xs:
   for _ in range(2):la.eigh(a,driver=d,check_finite=False)
   q=[]
   for _ in range(3):
    t=time.perf_counter();w,v=la.eigh(a,driver=d,check_finite=False);q.append(time.perf_counter()-t)
   ts.append(statistics.median(q));rr.append(np.max(abs(a@v-v*w)))
  out.append({"n":n,"driver":d,"median":statistics.median(ts),"total":sum(ts),"res":float(max(rr))})
open("/tmp/alt.json","w").write(json.dumps({"scipy":scipy.__version__,"results":out}))
