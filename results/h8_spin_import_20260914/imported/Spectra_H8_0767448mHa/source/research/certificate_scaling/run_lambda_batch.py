"""Bounded shared-instance experiment batch; no provider lifecycle actions."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/certificate_scaling/lambda_batch'
OUT.mkdir(parents=True,exist_ok=True)
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',CUDA_PATH='/usr')


def job(name,args,limit=180):
    start=time.monotonic()
    with (OUT/(name+'.log')).open('w') as log:
        try:
            proc=subprocess.run([sys.executable,*args],cwd=ROOT,env=ENV,
                                stdout=log,stderr=subprocess.STDOUT,timeout=limit)
            status={'returncode':proc.returncode,'status':'finished' if proc.returncode==0 else 'failed'}
        except subprocess.TimeoutExpired:
            status={'returncode':None,'status':'timeout'}
    status.update(name=name,command=args,timeout_seconds=limit,wall_seconds=time.monotonic()-start)
    (OUT/(name+'_status.json')).write_text(json.dumps(status,indent=2)+'\n')
    return status


if __name__=='__main__':
    started=time.monotonic()
    (ROOT/'results/certificate_scaling/operator_pricing').mkdir(parents=True,exist_ok=True)
    fixture='results/marginal_molecule_stress/h4_square_degree3_certificate.json'
    second='results/marginal_molecule_stress/accepted_degree3_certificate.json'
    prefix='results/certificate_scaling/lambda_batch/'
    receipts=[]
    # One CPU discovery search overlaps one GPU experiment; no competing GPU jobs.
    with ThreadPoolExecutor(max_workers=1) as pool:
        pricing=pool.submit(job,'pricing_m6',['-m','research.certificate_scaling.operator_pricing_sos','--modes','6','--budget','2'],300)
        for name,fx,backend in [('h4_gpu',fixture,'cupy'),('h4_cpu',fixture,'numpy'),('heldout_gpu',second,'cupy')]:
            receipts.append(job(name,['research/certificate_scaling/low_rank_compression.py','--campaign','--fixture',fx,'--outputdir',prefix+name,'--backend',backend]))
        receipts.append(job('locality_gpu',['research/certificate_scaling/locality_compression.py','--sites','8','16','32','64','--clusters','2','3','4','--couplings','0','1','2','4','8','--global-max-sites','0','--backend','cupy','--output',prefix+'locality_gpu.json']))
        receipts.append(job('sparse_factors',['-S','research/certificate_scaling/sparse_factor_campaign.py','--fixture',fixture,'--outputdir',prefix+'sparse_factors']))
        for rank in (8,32,64):
            receipts.append(job('ideal_repair_'+str(rank),['research/certificate_scaling/ideal_repair.py','--input',prefix+f'h4_gpu/rank_{rank}.json','--output',prefix+f'ideal_rank_{rank}.json','--receipt',prefix+f'ideal_rank_{rank}_receipt.json']))
        receipts.append(job('structured_replay',['-S','research/certificate_scaling/full_residual_replay.py',*[prefix+f'ideal_rank_{r}.json' for r in (8,32,64)],'--output',prefix+'structured_replay.json']))
        receipts.append(pricing.result())
    (OUT/'batch_status.json').write_text(json.dumps({'wall_seconds':time.monotonic()-started,'jobs':receipts},indent=2)+'\n')
    print(json.dumps({'jobs':len(receipts),'statuses':[(r['name'],r['status']) for r in receipts]}))
