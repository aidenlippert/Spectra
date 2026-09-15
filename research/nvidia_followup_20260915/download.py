"""Verify owned-host result bundles and import only new files or identical inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time
from research.nvidia_followup_20260915.remote import ROOT,OUT,KEY


def run(mode):
    start=time.monotonic()
    ready=json.loads((OUT/'ready.json').read_text())
    owned=json.loads((OUT/'instance.json').read_text())
    if ready['id']!=owned['id']:raise ValueError('Owned host changed')
    stem=mode+'_download'
    archive=OUT/(stem+'.tar.gz')
    if archive.exists():raise ValueError('Preserve existing downloaded archive')
    prefix=['scp','-i',str(KEY),'-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',
        '-o',f'UserKnownHostsFile={OUT}/known_hosts']
    origin=f"ubuntu@{ready['ip']}:/home/ubuntu/spectra-nvidia/"
    with (OUT/(stem+'.log')).open('x') as stream:
        p=subprocess.run(prefix+[origin+stem+'.tar.gz',origin+stem+'.sha256',str(OUT)],
            stdout=stream,stderr=subprocess.STDOUT,timeout=180)
    if p.returncode:raise RuntimeError('Download failed; inspect preserved log')
    digest=hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest!=(OUT/(stem+'.sha256')).read_text().strip():raise ValueError('Archive digest mismatch')
    dest=OUT/('download_'+mode)
    dest.mkdir(exist_ok=False)
    with tarfile.open(archive) as stream:stream.extractall(dest,filter='data')
    manifest=json.loads((dest/(stem+'_manifest.json')).read_text())
    imported=0
    for name,expected in manifest.items():
        if Path(name).is_absolute() or '..' in Path(name).parts:raise ValueError('Unsafe bundle entry')
        file=dest/name
        if hashlib.sha256(file.read_bytes()).hexdigest()!=expected:raise ValueError(('Downloaded file changed',name))
        if name.startswith(('results/transfer_solver_20260915/adaptive/h10_correlated_guide/',
            'results/transfer_solver_20260915/cases/h8_cold/')):
            target=ROOT/name
            if target.exists():
                if target.read_bytes()!=file.read_bytes():raise ValueError(('Refuse changed local input',name))
            else:
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(file,target)
                imported+=1
    receipt={'archive_sha256':digest,'files_verified':len(manifest),'new_files_imported':imported,
        'wall_seconds':time.monotonic()-start,'existing_inputs_changed':False}
    with (OUT/(stem+'_receipt.json')).open('x') as stream:json.dump(receipt,stream,indent=2)
    print(json.dumps(receipt))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('mode',choices=('proposal','final'))
    a=p.parse_args();run(a.mode)
