"""Give each proposal a separate original-model, standard-library replay folder."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from research.acceptance_channels_20260915.campaign import ROOT,OUT,CASE,dump
from research.interacting_scaling_20260915.complete import run as complete


def run(tag,case=CASE,name=None):
    proposal=case/tag
    witness=proposal/'export/certificate.json'
    if not witness.is_file():raise ValueError('A complete exported proposal is required')
    before=hashlib.sha256(witness.read_bytes()).hexdigest()
    replay=OUT/'replays'/(name or tag);replay.mkdir(parents=True,exist_ok=False)
    dependencies={}
    for name in ('fixture.json','upper.json','nonsinglet.json','rotation.json','mps'):
        source=(case/name).resolve(strict=True)
        if source.is_dir():
            shutil.copytree(source,replay/name)
        else:
            shutil.copyfile(source,replay/name)
            dependencies[name]=hashlib.sha256(source.read_bytes()).hexdigest()
    (replay/'proposal').symlink_to(proposal.resolve(strict=True),target_is_directory=True)
    dump(replay/'dependencies.json',{'source_proposal':str(proposal),'proposal_sha256':before,
        'source_sha256':dependencies,'prior_discovery_and_preparation_additional':True,
        'replay_is_not_a_cold_calculation':True})
    complete(ROOT/'results/interacting_scaling_20260915/models/h12_heldout',
        replay,'proposal',replay/'exact')
    if hashlib.sha256(witness.read_bytes()).hexdigest()!=before:
        raise ValueError('Proposal changed during replay')
    result=json.loads((replay/'original_interval.json').read_text())
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('tag');p.add_argument('--case',type=Path,default=CASE)
    p.add_argument('--name');a=p.parse_args()
    if Path(a.tag).name!=a.tag or a.tag in ('.','..'):p.error('One local proposal name required')
    if a.name and (Path(a.name).name!=a.name or a.name in ('.','..')):p.error('One new replay name required')
    run(a.tag,a.case.resolve(),a.name)
