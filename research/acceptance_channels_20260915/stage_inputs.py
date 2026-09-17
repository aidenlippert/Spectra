"""Materialize byte-identical read-only-use inputs outside the synced Documents tree."""
from concurrent.futures import ThreadPoolExecutor,as_completed
import hashlib
import json
from pathlib import Path
import time
from research.acceptance_channels_20260915.campaign import ROOT,OUT,CASE,dump


def main():
    started=time.monotonic()
    cache=Path('/Users/aidenlippert/.cache/spectra/acceptance_channels_20260915')
    cache.mkdir(parents=True,exist_ok=False)
    work=OUT/'h12_cached';work.mkdir(exist_ok=False)
    jobs=[]
    for source in sorted((CASE/'prepared').iterdir()):
        if source.is_file():jobs.append((source,cache/'prepared'/source.name))
    for name in ('fixture.json','upper.json','nonsinglet.json','rotation.json','design.json'):
        jobs.append((CASE/name,cache/name))
    for source in (CASE/'mps').rglob('*'):
        if source.is_file():jobs.append((source,cache/'mps'/source.relative_to(CASE/'mps')))
    for source,name in ((CASE/'dense_t2/checkpoint.npz','dense4_checkpoint.npz'),
            (OUT/'dense_channels.json','dense_channels.json'),(OUT/'dense_channels_8.json','dense_channels_8.json')):
        jobs.append((source,cache/name))
    def copy_one(job):
        source,target=job;target.parent.mkdir(parents=True,exist_ok=True)
        error=None
        for attempt in range(2):
            try:
                data=source.read_bytes();break
            except OSError as caught:error=caught
        else:raise error
        digest=hashlib.sha256(data).hexdigest()
        with target.open('xb') as stream:stream.write(data)
        if hashlib.sha256(target.read_bytes()).hexdigest()!=digest:
            raise ValueError('Input cache copy failed hash comparison')
        return {'source':str(source.resolve()),'cache':str(target),'bytes':len(data),'sha256':digest}
    receipts=[]
    with ThreadPoolExecutor(max_workers=8) as pool:
        for future in as_completed([pool.submit(copy_one,job) for job in jobs]):
            receipts.append(future.result())
            if len(receipts)%50==0:print(json.dumps({'copied':len(receipts),'total':len(jobs)}),flush=True)
    expected=json.loads((ROOT/'results/interacting_scaling_20260915/cases/h12_main_refined/continuation_dependency.json').read_text())['source_sha256']
    checked=0
    for record in receipts:
        p=Path(record['cache'])
        if p.parent.name=='prepared':
            match=[v for k,v in expected.items() if k.endswith('/prepared/'+p.name)]
            if len(match)!=1 or record['sha256']!=match[0]:raise ValueError('Cached prepared input differs from sealed dependency')
            checked+=1
    for name in ('prepared','fixture.json','upper.json','nonsinglet.json','rotation.json','design.json','mps'):
        source=cache/name;(work/name).symlink_to(source,target_is_directory=source.is_dir())
    dump(OUT/'input_cache.json',{'seconds':time.monotonic()-started,'cache':str(cache),'case':str(work),
        'copied_files':len(receipts),'copied_bytes':sum(r['bytes'] for r in receipts),
        'sealed_prepared_hashes_matched':checked,'files':receipts,
        'source_files_unchanged':True,'not_an_algorithmic_speedup':True})
    print(json.dumps({'case':str(work),'seconds':time.monotonic()-started,'files':len(receipts)}),flush=True)


if __name__=='__main__':main()
