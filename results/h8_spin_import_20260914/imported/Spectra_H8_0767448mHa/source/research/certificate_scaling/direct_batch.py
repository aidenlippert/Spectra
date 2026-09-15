"""Bounded independent process experiments for a warm two-host Lambda pool."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
FIXTURES={'square':'results/marginal_molecule_stress/h4_square_degree3_certificate.json',
          'rectangle':'results/marginal_molecule/h4_compact_certificate.json',
          'h6':'results/marginal_h6/fixture.json'}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--lane',choices=['lp','compare'],required=True)
    p.add_argument('--workers',type=int,default=3)
    p.add_argument('--outputdir',type=Path,required=True)
    a=p.parse_args()
    if not 1<=a.workers<=4: raise ValueError('one to four concurrent jobs per host')
    out=a.outputdir.resolve();out.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for name,fixture in FIXTURES.items():
        if a.lane=='lp':
            for body in (1,2):
                for budget in (0,128,512,2048,8192):
                    label=f'{name}_body{body}_a{budget}'
                    jobs.append((label,[sys.executable,'research/certificate_scaling/direct_sparse_discovery.py',
                        '--fixture',fixture,'--max-atoms',str(budget),'--ideal-body',str(body),
                        '--time-limit','90','--outputdir',str(out/label)]))
        else:
            jobs.append((f'{name}_quadratic',[sys.executable,'research/certificate_scaling/direct_quadratic_baseline.py',
                         '--fixture',fixture,'--outputdir',str(out/f'{name}_quadratic')]))
            for budget in (256,1024):
                label=f'{name}_pair{budget}'
                jobs.append((label,[sys.executable,'research/certificate_scaling/direct_pair_sdp.py',
                    '--fixture',fixture,'--budget',str(budget),'--outputdir',str(out/label)]))
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(job):
        label,cmd=job; start=time.monotonic()
        with (out/f'{label}.log').open('w') as f:
            try:
                r=subprocess.run(cmd,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
                status={'exit_code':r.returncode,'timed_out':False}
            except subprocess.TimeoutExpired:
                status={'exit_code':None,'timed_out':True}
        return dict(label=label,command=cmd,wall_seconds=time.monotonic()-start,**status)
    started=time.monotonic(); receipts=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            r=future.result();receipts.append(r)
            print(json.dumps(r),flush=True)
            (out/'jobs.json').write_text(json.dumps(receipts,indent=2)+'\n')
    manifest={str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest()
              for f in out.rglob('*') if f.is_file() and f.name!='manifest.json'}
    (out/'manifest.json').write_text(json.dumps({'files':manifest,'total_wall_seconds':time.monotonic()-started,
        'finished_utc':datetime.now(timezone.utc).isoformat(),'lane':a.lane,'workers':a.workers},indent=2)+'\n')

if __name__=='__main__':main()
