"""Threshold audit for conditioned H6 row-column certificate."""
import copy,json,time,hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
from experiments.marginal_symbolic import verify
SRC=Path('results/lambda_runs/cubic_precision/downloaded/precision_scs/h6_row-column/certificate.json'); SRC_RECEIPT=SRC.parent/'receipt.json'; OUT=Path('results/certificate_scaling/cubic_precision/factor_structure_conditioned')
def run():
 d=json.loads(SRC.read_text()); norms=[np.linalg.norm(r) for b in d['blocks'] for r in b['factor']]; trials=[]
 for th in (1e5,2e5,5e5,1e6,2e6):
  x=copy.deepcopy(d); kept=nz=0; loss=0
  for b in x['blocks']:
   rr=[]
   for r in b['factor']:
    n=float(np.linalg.norm(r))
    if n>=th: rr.append(r);kept+=1;nz+=sum(v!=0 for v in r)
    else: loss+=n*n
   b['factor']=rr
  p=OUT/f'threshold_{int(th)}.json'; OUT.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,separators=(',',':')))
  st=time.perf_counter()
  try: e=verify(x); status='verified_lower'
  except Exception as exc:e={'error':str(exc)};status='rejected'
  up=Fraction('-63330586263215870660998836114091/10000000000139879611500000000000'); lo=Fraction(e['lower']) if status=='verified_lower' else None; width=up-lo if lo is not None else None
  trials.append({'threshold':th,'kept_rows':kept,'nonzeros':nz,'bytes':p.stat().st_size,'discarded_unscaled_row_l2_squared':loss,'status':status,'exact':e,'upper_witness':str(up),'lower':str(lo) if lo is not None else None,'interval_width':str(width) if width is not None else None,'passes_accuracy_target':bool(width is not None and width<=Fraction(16,10000)),'wall_seconds':time.perf_counter()-st})
 eligible=[t for t in trials if t['passes_accuracy_target']]; best=min(eligible,key=lambda t:t['kept_rows'],default=None)
 sr=json.loads(SRC_RECEIPT.read_text()); receipt={'source':str(SRC),'source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'upper_witness':str(up),'original_discovery_wall_seconds':sr.get('wall_seconds'),'trials':trials,'best':best,'selection_rule':'fewest retained rows among exact replays with exact Fraction interval width <= 16/10000'}
 (OUT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');
 if best: (OUT/'compressed_certificate.json').write_bytes((OUT/f"threshold_{int(best['threshold'])}.json").read_bytes())
 return receipt
if __name__=='__main__': print(json.dumps(run(),indent=2))
