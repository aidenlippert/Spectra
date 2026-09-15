"""Three-threshold H8 post-discovery factor compression plus wedge replay."""
import copy,json,time,hashlib,subprocess,sys
from pathlib import Path
import numpy as np
from fractions import Fraction as F
SRC=Path('results/lambda_runs/cubic_precision/downloaded/h8_scs_quotient_conditioned/certificate.json');OUT=Path('results/certificate_scaling/cubic_precision/factor_structure_h8');UP=F(json.loads(Path('results/certificate_scaling/cubic_precision/reference_h8_full_checked/h8/upper.json').read_text())['upper'])
def main():
 d=json.loads(SRC.read_text()); trials=[]
 for th in (1e5,5e5,1e6):
  x=copy.deepcopy(d); rows=nz=0
  for b in x['blocks']:
   keep=[]
   for r in b['factor']:
    if np.linalg.norm(r)>=th: keep.append(r);rows+=1;nz+=sum(v!=0 for v in r)
   b['factor']=keep
  p=OUT/f'threshold_{int(th)}.json';OUT.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,separators=(',',':')))
  o=OUT/f'threshold_{int(th)}_wedge';st=time.perf_counter();r=subprocess.run([sys.executable,'-m','research.certificate_scaling.wedge_spectral_bound','--certificate',str(p),'--out',str(o)],capture_output=True,text=True);wall=time.perf_counter()-st
  if r.returncode:raise RuntimeError(r.stderr[-1000:])
  rec=json.loads((o/'receipt.json').read_text())
  if rec['certificate_sha256']!=hashlib.sha256(p.read_bytes()).hexdigest():raise ValueError('Stale trial receipt')
  proof=json.loads((o/'witness.json').read_text());lo=F(rec['lower']);width=UP-lo
  trials.append({'threshold':th,'rows':rows,'nonzeros':nz,'certificate_bytes':p.stat().st_size,'spectral_receipt':rec,'spectral_wall_seconds':wall,'lower':str(lo),'upper':str(UP),'interval_width':str(width),'passes_accuracy_target':0<=width<=F(16,10000),
                 'spectral_factor_rows':sum(len(c['factor'][0]) for b in proof['bodies'] for c in b['components']),
                 'spectral_proof_bytes':(o/'witness.json').stat().st_size})
 eligible=[t for t in trials if t['passes_accuracy_target']];best=min(eligible,key=lambda t:t['rows'],default=None)
 receipt={'source':str(SRC),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'original_discovery_wall_seconds':721.5632876,'upper':str(UP),'trials':trials,'best':best,'selection_rule':'fewest factor rows among wedge-replayed intervals <= 0.0016'}
 if best:
  (OUT/'compressed_certificate.json').write_bytes((OUT/f"threshold_{int(best['threshold'])}.json").read_bytes())
  cmd=[sys.executable,'-S','-m','research.certificate_scaling.cubic_interval_replay','--certificate',str(OUT/'compressed_certificate.json'),'--reference','results/certificate_scaling/cubic_precision/reference_h8_full_checked/h8/upper.json','--method','spectral','--proof',str(OUT/f"threshold_{int(best['threshold'])}_wedge"/'witness.json'),'--out',str(OUT/'independent_interval.json')]
  subprocess.run(cmd,check=True,capture_output=True,text=True)
  receipt['independent_interval']=json.loads((OUT/'independent_interval.json').read_text())
 receipt['all_trial_wall_seconds']=sum(t['spectral_wall_seconds'] for t in trials)
 (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
