"""Finalize accounting after all experiment runners finish, then seal bytes.

This is file inventory only, outside the numerical experiment wall ledger.
Future numerical work must use new output directories rather than this seal.
"""
import json
from pathlib import Path
from research.reconstruction_compression_20260914.inputs import ROOT,OUT,sha,dump
from research.reconstruction_compression_20260914.report import main

def seal():
    target=OUT/'manifest.json'
    if target.exists():raise ValueError('Campaign already sealed')
    main()
    files={}
    for folder in (ROOT/'research/reconstruction_compression_20260914',OUT):
        for p in sorted(folder.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p!=target:
                files[str(p.relative_to(ROOT))]={'bytes':p.stat().st_size,'sha256':sha(p)}
    dump(target,{'kind':'reconstruction_preserving_campaign_inventory_v1','files':files,
                 'scope':'Source and result files; excludes interpreter caches and this manifest itself.'})
    for p,v in files.items():
        if sha(ROOT/p)!=v['sha256']:raise ValueError('File changed during sealing')
    print(json.dumps({'sealed_files':len(files),'manifest_sha256':sha(target),'sealed_bytes':sum(v['bytes'] for v in files.values())},indent=2))

if __name__=='__main__':seal()
