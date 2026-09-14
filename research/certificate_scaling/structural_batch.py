"""Bounded canonical/localized scaling sweep on existing warm hosts."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import argparse,subprocess,os,sys,time,json,hashlib
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser();p.add_argument('--lane',choices=['canonical','localized'],required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);jobs=[]
 for n in [4,6,8,10]:
  folder='active_space_ladder' if a.lane=='canonical' else ('active_space_ladder_boys' if n<=6 else 'active_space_ladder_boys_large')
  fixture=f'results/certificate_scaling/{folder}/h{n}/fixture.json'
  for mode in ['full','adaptive']:
   label=f'h{n}_{mode}';cmd=[sys.executable,'research/certificate_scaling/adaptive_block_discovery.py','--fixture',fixture,'--outputdir',str(a.out/label),'--symmetry','--width','4','--rounds','5','--seconds','180']
   if mode=='full':cmd+=['--full']
   jobs.append((label,cmd,240,fixture))
 if a.lane=='canonical':
  for n in [4,6]:
   fixture=f'results/certificate_scaling/active_space_ladder/h{n}/fixture.json';label=f'h{n}_full_nosym'
   jobs.append((label,[sys.executable,'research/certificate_scaling/adaptive_block_discovery.py','--fixture',fixture,'--outputdir',str(a.out/label),'--full','--seconds','120'],180,fixture))
  for name,fixture in [('square','results/marginal_molecule_stress/h4_square_degree3_certificate.json'),('rectangle','results/marginal_molecule/h4_compact_certificate.json')]:
   label=name+'_mixed_full';jobs.append((label,[sys.executable,'research/certificate_scaling/adaptive_block_discovery.py','--fixture',fixture,'--outputdir',str(a.out/label),'--symmetry','--full','--family','mixed','--ideal-body','2','--seconds','180'],240,fixture))
 def execute(job):
  label,cmd,limit,fixture=job;start=time.monotonic()
  with (a.out/(label+'.log')).open('w') as f:
   try:r=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),timeout=limit);code=r.returncode
   except subprocess.TimeoutExpired:code='timeout'
  return {'label':label,'command':cmd,'fixture':fixture,'exit_code':code,'wall_seconds':time.monotonic()-start}
 rows=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for fut in as_completed([pool.submit(execute,j) for j in jobs]):
   row=fut.result();rows.append(row);print(json.dumps(row),flush=True);(a.out/'jobs.json').write_text(json.dumps(rows,indent=2)+'\n')
 manifest={str(f.relative_to(a.out)):hashlib.sha256(f.read_bytes()).hexdigest() for f in a.out.rglob('*') if f.is_file() and f.name!='manifest.json'}
 (a.out/'manifest.json').write_text(json.dumps({'files':manifest},indent=2)+'\n')
if __name__=='__main__':main()
