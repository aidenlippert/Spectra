"""Seal completed evidence; this does not mark the research objective solved."""
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/collective_completion_20260914'
SRC=ROOT/'research/collective_completion_20260914'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def run():
    if (OUT/'manifest.json').exists():raise RuntimeError('Already sealed')
    audit=json.loads((OUT/'audit.json').read_text());result=json.loads((OUT/'final_result.json').read_text())
    if audit['changed'] or audit['pending'] or not result['acceptance_tests']['passed']:raise ValueError('Unfinished evidence gates')
    for p in (OUT/'runs').glob('*.json'):
        if json.loads(p.read_text())['status']=='starting':raise ValueError('A process is still pending')
    files={}
    for root in (SRC,OUT):
        for p in sorted(root.rglob('*')):
            if not p.is_file() or '__pycache__' in p.parts or p in (OUT/'manifest.json',OUT/'seal.json'):continue
            files[str(p.relative_to(ROOT))]={'sha256':digest(p),'bytes':p.stat().st_size}
    parent=ROOT/'results/reconstruction_compression_20260914/manifest.json'
    manifest={'kind':'sealed_research_evidence_v1','created_UTC':datetime.now(timezone.utc).isoformat(),
        'parent_manifest':str(parent.relative_to(ROOT)),'parent_sha256':digest(parent),
        'inherited_files_verified':audit['preservation_checked'],'files':files,
        'file_count':len(files),'total_bytes':sum(v['bytes'] for v in files.values()),
        'objective_status':{'original_H6_compact_target':result['best_verified']['h6']['accepted_target'],
                            'fresh_H6_compact_target':result['best_verified']['fresh_h6']['accepted_target'],
                            'H8_compact_target':result['best_verified']['h8']['accepted_target'],
                            'general_many_body_solution':False}}
    path=OUT/'manifest.json'
    with path.open('x') as f:json.dump(manifest,f,indent=2,sort_keys=True);f.write('\n')
    receipt={'manifest_sha256':digest(path),'file_count':len(files),'total_bytes':manifest['total_bytes'],
             'status':'Evidence sealed; unresolved research goals remain explicitly recorded'}
    with (OUT/'seal.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':run()
