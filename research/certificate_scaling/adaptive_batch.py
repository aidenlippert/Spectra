"""Independent bounded jobs on the existing warm Lambda hosts."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys,time

ROOT=Path(__file__).resolve().parents[2]
FIX={'square':'results/marginal_molecule_stress/h4_square_degree3_certificate.json',
     'rectangle':'results/marginal_molecule/h4_compact_certificate.json',
     'h6':'results/marginal_h6/fixture.json'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--lane',choices=['quadratic','cubic','optimized_quadratic','optimized_cubic'],required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=3)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True);jobs=[]
    specs=[]
    if a.lane=='optimized_quadratic':
        specs=[(name,'quadratic',4,1,200,300) for name in FIX]
    elif a.lane=='optimized_cubic':
        specs=[(name,'mixed',16,2,100,300) for name in ('square','rectangle')]
        specs += [('h6','local',8,2,100,300)]
    elif a.lane=='quadratic':
        specs=[(name,'quadratic',w,1,80,180) for name in FIX for w in (2,4,8,0)]
    else:
        specs=[(name,'mixed',w,2,50,180) for name in ('square','rectangle') for w in (4,8,16,0)]
        specs += [('h6','local',w,2,30,240) for w in (8,0)]
    for name,family,width,body,its,seconds in specs:
        label=f'{name}_{family}_w{width}'
        cmd=[sys.executable,'research/certificate_scaling/adaptive_factor_pricing.py','--fixture',FIX[name],
             '--family',family,'--width',str(width),'--ideal-body',str(body),'--iterations',str(its),
             '--batch','64','--atom-cap','6500','--seconds',str(seconds),'--outputdir',str(out/label)]
        if a.lane.startswith('optimized_'):cmd += ['--half-rows','--prune']
        jobs.append((label,cmd,seconds+60))
    if a.lane=='quadratic':
        for name in ('square','rectangle'):
            label=name+'_full_pair'
            jobs.append((label,[sys.executable,'research/certificate_scaling/full_pair_cone.py',
                         '--fixture',FIX[name],'--outputdir',str(out/label)],150))
    def execute(job):
        label,cmd,timeout=job;start=time.monotonic()
        with (out/(label+'.log')).open('w') as f:
            try:r=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,
                        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),timeout=timeout);code=r.returncode
            except subprocess.TimeoutExpired:code='timeout'
        return {'label':label,'command':cmd,'exit_code':code,'wall_seconds':time.monotonic()-start}
    start=time.monotonic();records=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(execute,j) for j in jobs]):
            r=f.result();records.append(r);print(json.dumps(r),flush=True)
            (out/'jobs.json').write_text(json.dumps(records,indent=2)+'\n')
    manifest={str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest()
              for f in out.rglob('*') if f.is_file() and f.name!='manifest.json'}
    (out/'manifest.json').write_text(json.dumps({'files':manifest,'batch_seconds':time.monotonic()-start},indent=2)+'\n')

if __name__=='__main__':main()
