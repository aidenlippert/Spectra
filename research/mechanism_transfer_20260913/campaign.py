"""Sequential bounded template controls on frozen training and held-out cases."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/mechanism_transfer_20260913'
ENV={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1'}


def save(path,obj): path.write_text(json.dumps(obj,indent=2)+'\n')


def case(name,out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    from research.mechanism_transfer_20260913.discovery import Model
    from research.mechanism_transfer_20260913.rule import generate
    from research.molecular_collective_20260913.core import extract,propose_tail,tail_replay
    from research.certificate_scaling.streaming_reference_upper import upper
    choice_path=OUT/'frozen_rule.json';choice=json.loads(choice_path.read_text())
    if hashlib.sha256((Path(__file__).parent/'rule.py').read_bytes()).hexdigest()!=choice['rule_source_sha256']:
        raise ValueError('Frozen generation rule changed')
    if name=='h6_train':
        prior=ROOT/'results/molecular_collective_20260913/campaign/h6'
        fixture=prior/'fixture.json';tail_path=prior/'rank_10/tail.json'
        ref=ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h6/upper.json'
    elif name=='h8_heldout':
        prior=ROOT/'results/molecular_collective_20260913/campaign/h8'
        fixture=prior/'fixture.json';tail_path=prior/'rank_14/tail.json'
        ref=ROOT/'results/certificate_scaling/active_space_ladder_references_aligned/h8/upper.json'
    else:
        prior=OUT/'fixtures'/name;fixture=prior/'fixture.json';ref=prior/'upper.json';tail_path=None
    data=json.loads(fixture.read_text());reference=json.loads(ref.read_text())
    U,reference_receipt=upper(data,reference['independent_upper'])
    if U!=F(reference['upper']): raise ValueError('Upper reference did not replay')
    before=time.monotonic()
    if tail_path is not None:
        tail=json.loads(tail_path.read_text());tail_discovery={'scope':'Frozen existing tail; prior discovery remains in its ledger.'}
    else:
        p=extract(data);tail,tail_discovery=propose_tail(data,p,2*p['spatial']-2)
    tail_receipt=tail_replay(data,tail)
    save(out/'fixture.json',data);save(out/'tail.json',tail);save(out/'upper.json',reference)
    metadata={'case':name,'upper_Ha':str(U),'reference_replay':reference_receipt,
        'tail':tail_receipt,'tail_discovery':tail_discovery,'tail_preparation_seconds':time.monotonic()-before,
        'frozen_rule_sha256':hashlib.sha256(choice_path.read_bytes()).hexdigest(),'budget_seconds':240.,
        'source_fixture_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest(),'trials':[]}
    save(out/'summary.json',metadata)
    model=Model(data,tail);span,rule=generate(model,tuple(choice['powers']))
    metadata.update(model_construction=model.construction,rule_construction=rule)
    save(out/'generated_span.json',span);save(out/'summary.json',metadata)
    for arm,selected,separate in (('quadratic',[],False),('coupled',span,False),('separate',span,True)):
        if time.monotonic()-start>240-35:
            metadata['unrun_arm']=arm;break
        trial_start=time.monotonic();metadata['pending_arm']=arm;save(out/'summary.json',metadata)
        try:
            receipt,_=model.solve(selected,out/arm,240-(time.monotonic()-start),separate=separate,solver_cap=45.)
            L=F(receipt['accepted']['original_lower_Ha']);done=time.monotonic()-start
            record={'arm':arm,'accepted':True,'complete_at_seconds':done,'inside_budget':done<=240,
                'lower_Ha':str(L),'width_Ha':str(U-L),'width_mHa':float(1000*(U-L)),
                'trial_seconds':time.monotonic()-trial_start,'directions':receipt['directions']}
        except Exception as error:
            record={'arm':arm,'accepted':False,'error':str(error),'trial_seconds':time.monotonic()-trial_start}
        metadata['trials'].append(record);metadata['pending_arm']=None;save(out/'summary.json',metadata)
        print(json.dumps(record),flush=True)
    eligible=[r for r in metadata['trials'] if r['accepted'] and r['inside_budget']]
    metadata.update(best_arm=min(eligible,key=lambda r:F(r['width_Ha']))['arm'] if eligible else None,
        wall_seconds=time.monotonic()-start,complete=True)
    save(out/'summary.json',metadata)


def run():
    OUT.mkdir(exist_ok=True);records=[]
    for name,spacing in (('h6_1p2','1.2'),('h6_1p8','1.8')):
        dest=OUT/'fixtures'/name;dest.parent.mkdir(exist_ok=True);log=dest.with_suffix('.log')
        if dest.exists() or log.exists(): raise FileExistsError('Existing fixture work must be preserved')
        start=time.monotonic();failure=None
        with log.open('w') as stream:
            try:
                p=subprocess.run([sys.executable,'-m','research.mechanism_transfer_20260913.fixtures',
                    '--spacing',spacing,'--out',str(dest)],cwd=ROOT,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=120)
                if p.returncode: failure=f'exit_{p.returncode}'
            except subprocess.TimeoutExpired: failure='wall_timeout'
        row={'case':name,'phase':'fixture_generation','budget_seconds':120.,'wall_seconds':time.monotonic()-start,'failure':failure}
        records.append(row);save(OUT/'watchdogs.json',records);print(json.dumps(row),flush=True)
    for name in ('h6_train','h6_1p2','h6_1p8','h8_heldout'):
        dest=OUT/'campaign'/name;dest.parent.mkdir(exist_ok=True);log=dest.with_suffix('.log')
        if dest.exists() or log.exists(): raise FileExistsError('Existing case work must be preserved')
        start=time.monotonic();failure=None
        with log.open('w') as stream:
            try:
                p=subprocess.run([sys.executable,'-m','research.mechanism_transfer_20260913.campaign',
                    '--case',name,'--out',str(dest)],cwd=ROOT,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=240)
                if p.returncode:failure=f'exit_{p.returncode}'
            except subprocess.TimeoutExpired:failure='wall_timeout'
        row={'case':name,'phase':'discovery','budget_seconds':240.,'wall_seconds':time.monotonic()-start,'failure':failure}
        records.append(row);save(OUT/'watchdogs.json',records);print(json.dumps(row),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case');ap.add_argument('--out',type=Path);args=ap.parse_args()
    if args.case:case(args.case,args.out)
    else:run()
