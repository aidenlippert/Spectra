"""Bounded size/dispersion/frontier campaign with all discovery and replay costs."""
from pathlib import Path
import argparse
import json
import time

from research.collective_interference_20260913.collective import model,propose,replay


def run(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);cases=[]
    for length in (1,2,4,8,16,32,64):cases.append((f'exact_L{length}',model(length,spread='0'),1))
    for length in (2,4,8,16,32,64):cases.append((f'dispersed_L{length}',model(length),1))
    for length in (2,4,8,16):cases.append((f'refined_L{length}',model(length),2))
    for gap in ('1','1/4','1/20'):cases.append((f'gap_{gap.replace("/","_")}',model(8,gap=gap),1))
    for coupling in ('0','6/5'):cases.append((f'coupling_{coupling.replace("/","_")}',model(16,hybridization=coupling),1))
    rows=[]
    for name,data,bins in cases:
        t=time.monotonic();directory=out/name;directory.mkdir()
        certificate,discovery=propose(data,bins);accepted=replay(data,certificate)
        (directory/'model.json').write_text(json.dumps(data,indent=2)+'\n')
        raw=json.dumps(certificate,separators=(',',':'))+'\n';(directory/'certificate.json').write_text(raw)
        row={'name':name,'bins_per_band':bins,'model_bytes':len(json.dumps(data,separators=(',',':')).encode()),
             'certificate_bytes':len(raw.encode()),'discovery':discovery,'accepted':accepted,
             'case_wall_seconds':time.monotonic()-t}
        (directory/'receipt.json').write_text(json.dumps(row,indent=2)+'\n');rows.append(row)
        print(json.dumps({'name':name,'physical_modes':accepted['physical_modes'],'retained_modes':accepted['retained_modes'],
                          'width':accepted['width_float'],'seconds':row['case_wall_seconds']}),flush=True)
    result={'rows':rows,'wall_seconds':time.monotonic()-start,
            'discovery_seconds':sum(x['discovery']['proposal_seconds'] for x in rows),
            'accepting_replay_seconds':sum(x['accepted']['replay_seconds'] for x in rows),
            'scope':'Structured family, abstract energy units. No physical many-body sector enumerated in discovery or replay.'}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.out)
