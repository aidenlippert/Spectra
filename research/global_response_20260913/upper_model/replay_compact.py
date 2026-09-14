import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
from research.global_response_20260913.slater import check
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'results/global_response_20260913/upper_model'
rows=json.loads((OUT/'compact_certificates.json').read_text()); out=[]
for x in rows:
 d=json.loads((ROOT/x['fixture']).read_text()); c=x['certificate']; out.append({'case':x['case'],'replay':check(d,c)})
(OUT/'compact_replay.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
