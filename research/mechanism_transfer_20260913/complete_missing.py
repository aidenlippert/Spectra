"""Recover only geometry cases that stopped before producing inputs."""
from pathlib import Path
import json
import subprocess
import sys
import time
from research.mechanism_transfer_20260913.campaign import OUT,ROOT,ENV,save


def run():
    original=json.loads((OUT/'watchdogs.json').read_text());records=[]
    for name,spacing in (('h6_1p2','1.2'),('h6_1p8','1.8')):
        dest=OUT/'fixtures'/name
        if dest.exists():raise FileExistsError('Recovery only admits missing geometry inputs')
        for phase in ('fixture_generation','discovery'):
            spent=sum(r['wall_seconds'] for r in original if r['case']==name and r['phase']==phase)
            limit=(120 if phase=='fixture_generation' else 240)-spent
            if phase=='fixture_generation':
                command=[str(ROOT/'.venv-molecule/bin/python'),'-m','research.mechanism_transfer_20260913.fixtures',
                    '--spacing',spacing,'--out',str(dest)]
                log=dest.with_name(name+'_runtime_fixed.log')
            else:
                target=OUT/'campaign'/(name+'_inputs_fixed')
                if target.exists():raise FileExistsError('Recovery output already exists')
                command=[sys.executable,'-m','research.mechanism_transfer_20260913.campaign',
                    '--case',name,'--out',str(target)]
                log=target.with_suffix('.log')
            if log.exists():raise FileExistsError('Preserve the recovery log')
            start=time.monotonic();failure=None
            with log.open('w') as stream:
                try:
                    proc=subprocess.run(command,cwd=ROOT,env=ENV,stdout=stream,stderr=subprocess.STDOUT,timeout=limit)
                    if proc.returncode:failure=f'exit_{proc.returncode}'
                except subprocess.TimeoutExpired:failure='wall_timeout'
            row={'case':name,'phase':phase,'remaining_budget_seconds':limit,
                'previous_failed_seconds':spent,'wall_seconds':time.monotonic()-start,'failure':failure,
                'reason':'The original fixture attempt used the numerical-SDP Python without PySCF; use the existing molecule runtime. No scientific result or rule was changed.'}
            records.append(row);save(OUT/'recovery_watchdogs.json',records);print(json.dumps(row),flush=True)


if __name__=='__main__':run()
