"""Byte-verified local cache of the installed SciPy package; no version change."""
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
from pathlib import Path
import time
from research.acceptance_channels_20260915.campaign import ROOT,OUT,dump


def main():
    started=time.monotonic()
    source=ROOT/'.venv-correlated/lib/python3.12/site-packages'
    destination=Path('/Users/aidenlippert/.cache/spectra/acceptance_channels_20260915/runtime')
    destination.mkdir(parents=True,exist_ok=False)
    paths=[p for folder in (source/'scipy',source/'scipy-1.18.1.dist-info') for p in folder.rglob('*') if p.is_file()]
    def copy(path):
        data=path.read_bytes();target=destination/path.relative_to(source)
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        digest=hashlib.sha256(data).hexdigest()
        if hashlib.sha256(target.read_bytes()).hexdigest()!=digest:raise ValueError('Runtime cache hash mismatch')
        return {'file':str(path.relative_to(source)),'bytes':len(data),'sha256':digest}
    records=[]
    with ThreadPoolExecutor(max_workers=8) as pool:
        for future in as_completed([pool.submit(copy,p) for p in paths]):
            records.append(future.result())
            if len(records)%250==0:print(f'{len(records)}/{len(paths)} files copied',flush=True)
    dump(OUT/'runtime_cache.json',{'seconds':time.monotonic()-started,'version':'1.18.1','source':str(source),
        'destination':str(destination),'files':records,'installed_files_unchanged':True,'algorithmic_speedup_claimed':False})
    print(f'Cached {len(records)} installed files in {time.monotonic()-started:.3f} seconds',flush=True)


if __name__=='__main__':main()
