"""Stdlib exact interval replay; references supplied only after discovery."""
from pathlib import Path
from fractions import Fraction as F
import json,time,sys,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from experiments.marginal_symbolic import verify
from experiments.marginal_determinant_tree import DeterminantOracle
from research.certificate_scaling.direct_results_replay import witness

def replay(path,out,label,localized=False):
 cert=json.loads(path.read_text());name=path.parent.name.split('_')[0]
 if name in ['square','rectangle']:w=witness(name);reference='previous independent H4 reference';reference_seconds=None
 else:
  folder='active_space_ladder_references_aligned_localized' if localized else 'active_space_ladder_references_aligned'
  reference=Path('results/certificate_scaling',folder,name,'upper.json');w=json.loads(reference.read_text())['independent_upper']
  r=json.loads(reference.with_name('receipt.json').read_text());reference_seconds=r['elapsed_seconds']
 cert['sparse_independent_upper']=w;start=time.monotonic();low=F(verify(cert)['lower']);up=DeterminantOracle(cert).upper(w);assert up>=low
 fp=out/(label+'.json');fp.write_text(json.dumps(cert,separators=(',',':'))+'\n')
 return {'label':label,'lower':str(low),'upper':str(up),'width':str(up-low),'lower_float':float(low),'upper_float':float(up),
    'width_float':float(up-low),'passes_0_0016_Ha':up-low<=F(1,625),'reference':str(reference),'reference_generation_seconds':reference_seconds,
    'reference_support':len(w['states']),'replay_seconds':time.monotonic()-start,'certificate':str(fp),'certificate_bytes':fp.stat().st_size,
    'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'certificate_sha256':hashlib.sha256(fp.read_bytes()).hexdigest()}
if __name__=='__main__':
 root=Path('results/lambda_runs/structural_scaling/downloaded');out=Path('results/certificate_scaling/structural_intervals');out.mkdir(parents=True,exist_ok=True);rows=[]
 for lane in ['canonical','localized','cubic_small','cubic_large']:
  for path in sorted((root/lane).rglob('certificate.json')):rows.append(replay(path,out,lane+'_'+path.parent.name,lane=='localized'))
 repair=Path('results/certificate_scaling/structural_repair/h6_mixed_full/certificate.json')
 if repair.exists():rows.append(replay(repair,out,'repaired_h6_mixed_full'))
 (out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps({'replayed':len(rows),'passing':sum(r['passes_0_0016_Ha'] for r in rows),'seconds':sum(r['replay_seconds'] for r in rows)}))
