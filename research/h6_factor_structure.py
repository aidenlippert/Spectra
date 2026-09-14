"""Post-discovery factor-row compression audit for the downloaded H6 proof."""
import copy,json,time
from pathlib import Path
import numpy as np
from experiments.marginal_symbolic import verify

SRC=Path('results/lambda_runs/cubic_precision/downloaded/precision_scs/h6_none/certificate.json')
OUT=Path('results/certificate_scaling/cubic_precision/factor_structure')

def audit():
    cert=json.loads(SRC.read_text()); rows=[]; dimensions=[]
    for bi,b in enumerate(cert['blocks']):
        a=np.asarray(b['factor'],dtype=float); dimensions.append(len(a))
        for ri,r in enumerate(a): rows.append((float(np.linalg.norm(r)),bi,ri,float(np.sum(np.abs(r))),int(np.count_nonzero(r))))
    rows.sort(); total=sum(x[0]**2 for x in rows)
    receipt={'source':str(SRC),'blocks':len(cert['blocks']),'factor_rows':len(rows),'factor_nonzeros':sum(x[4] for x in rows),'dimensions':dimensions,'norm_l2_sum':total,'row_norm_quantiles':{str(q):float(np.quantile([x[0] for x in rows],q)) for q in (0,.01,.1,.5,.9,.99,1)},'tests':[]}
    for threshold in (0.0,1e4,1e5,1e6):
        trial=copy.deepcopy(cert); kept=0; nz=0; loss=0.
        for b in trial['blocks']:
            new=[]
            for r in b['factor']:
                n=float(np.linalg.norm(r));
                if n>=threshold or threshold==0: new.append(r); kept+=1; nz+=sum(v!=0 for v in r)
                else: loss+=n*n
            b['factor']=new
        started=time.perf_counter();
        try: exact=verify(trial); status='accepted'
        except Exception as e: exact={'error':str(e)}; status='rejected'
        receipt['tests'].append({'threshold':threshold,'kept_rows':kept,'nonzeros':nz,'discarded_unscaled_row_l2_squared':loss,'bytes':len(json.dumps(trial).encode()),'wall_seconds':time.perf_counter()-started,'status':status,'exact':exact})
    OUT.mkdir(parents=True,exist_ok=True); (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); return receipt

def write_best():
    cert=json.loads(SRC.read_text()); threshold=1e6; kept=nz=0; loss=0.0
    for b in cert['blocks']:
        rows=[]
        for r in b['factor']:
            n=float(np.linalg.norm(r))
            if n>=threshold: rows.append(r); kept+=1; nz+=sum(v!=0 for v in r)
            else: loss+=n*n
        b['factor']=rows
    OUT.mkdir(parents=True,exist_ok=True); path=OUT/'compressed_certificate.json'; path.write_text(json.dumps(cert,separators=(',',':'))+'\n')
    started=time.perf_counter(); exact=verify(cert); elapsed=time.perf_counter()-started
    lower=float(exact['lower_float']); upper=-6.333058626233; width=upper-lower
    receipt={'status':'verified_lower','source':str(SRC),'compressed_certificate':str(path),'threshold':threshold,'kept_rows':kept,'nonzeros':nz,'bytes':path.stat().st_size,'discarded_unscaled_row_l2_squared':loss,'exact_replay_seconds':elapsed,'exact_replay':exact,'upper_witness':upper,'interval_width':width,'passes_accuracy_target':width<=0.0016,'post_discovery_source_wall_seconds':259.6397602069992,'compression_and_replay_wall_seconds':elapsed}
    (OUT/'compressed_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n'); return receipt

if __name__=='__main__':
    import sys
    print(json.dumps(write_best() if '--best' in sys.argv else audit(),indent=2))
