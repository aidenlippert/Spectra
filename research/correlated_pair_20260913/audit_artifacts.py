"""Preservation and measured-cost inventory; no solver or numerical imports."""
import hashlib
import json
from pathlib import Path
import time


def run():
    start=time.monotonic();base=Path('results/correlated_pair_20260913');prior=json.loads((base/'preservation_before.json').read_text());changes=[]
    for name,expected in prior['files'].items():
        path=Path(name)
        if not path.exists():changes.append({'path':name,'missing':True});continue
        raw=path.read_bytes();current={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
        if current!=expected:changes.append({'path':name,'expected':expected,'current':current})
    receipts=[];unavailable=[]
    for path in sorted((base/'runs').glob('*.json')):
        try:row=json.loads(path.read_text())
        except ValueError:unavailable.append({'path':str(path),'bytes':path.stat().st_size,'reason':'Empty early optimizer output; timing receipt unavailable.'});continue
        if 'status' not in row:
            unavailable.append({'path':str(path),'reason':'Early optimizer summary without a durable runner receipt; partial internal timing only.','partial_internal_seconds':row.get('seconds')});continue
        if row.get('name','').startswith('preservation_and_accounting'):continue # Runner finishes after this script exits.
        receipts.append(row)
    result={'preserved_file_count':len(prior['files']),'changes_since_campaign_start':changes,
            'preexisting_parent_mismatches':prior.get('preexisting_parent_mismatches',[]),
            'logged_completed_process_wall_seconds':sum(r.get('wall_seconds',0) for r in receipts),
            'peak_logged_child_RSS_bytes':max(r.get('peak_child_RSS_bytes',0) for r in receipts),
            'counts_by_status':{s:sum(r['status']==s for r in receipts) for s in sorted({r['status'] for r in receipts})},
            'pending_experiments':[r['name'] for r in receipts if r['status']=='starting'],
            'missing_attempt_accounting':unavailable,'cost_scope':'Sum of logged process wall times, including failed and timed-out runs; not a complete end-to-end total because two early failed attempt timings are unavailable.',
            'historical_cubic_discovery_seconds':{'h6':259.309,'h8':721.563},
            'historical_cubic_cost_scope':'Retained historical measurements, excluded from the new campaign total; unsuccessful historical attempts are not fully accounted.',
            'preservation_scan_seconds':time.monotonic()-start}
    (base/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    if changes or result['pending_experiments']:raise ValueError('Preservation change or unfinished experiment')


if __name__=='__main__':run()
