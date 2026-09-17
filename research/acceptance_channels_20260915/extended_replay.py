"""Explicitly extended time envelope for the unchanged complete exact checker."""
import argparse
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time
from research.acceptance_channels_20260915.campaign import ROOT,OUT


def run(name,tag,case,output_name):
    for value in (name,tag,output_name):
        if not value or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for c in value):
            raise ValueError('Safe unique run and proposal names required')
    command=['/opt/homebrew/Caskroom/miniconda/base/bin/python','-B','-S','-m',
        'research.acceptance_channels_20260915.replay_candidate',tag,'--case',str(case),'--name',output_name]
    record={'name':name,'command':command,'timeout_seconds':1800,'status':'starting',
        'started_UTC':datetime.now(timezone.utc).isoformat(),
        'scope':'Same original complete exact checker; declared extended replay budget after the 900-second timeout.'}
    receipt=OUT/'runs'/(name+'.json');log=OUT/'runs'/(name+'.log')
    with receipt.open('x') as stream:json.dump(record,stream,indent=2)
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS'):env[key]='1'
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with log.open('x') as stream:
        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        record['pid']=process.pid
        receipt.write_text(json.dumps(record,indent=2)+'\n')
        try:
            code=process.wait(timeout=1800);record['status']='passed' if code==0 else 'failed'
        except subprocess.TimeoutExpired:
            record['status']='timeout';os.killpg(process.pid,signal.SIGTERM)
            try:code=process.wait(timeout=3)
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);code=process.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    record.update(exit_code=code,wall_seconds=time.monotonic()-start,
        child_user_seconds=after.ru_utime-before.ru_utime,child_system_seconds=after.ru_stime-before.ru_stime,
        peak_child_RSS_bytes=after.ru_maxrss*(1 if sys.platform=='darwin' else 1024),
        finished_UTC=datetime.now(timezone.utc).isoformat())
    receipt.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
    return 0 if record['status']=='passed' else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('tag');p.add_argument('case',type=Path)
    p.add_argument('output_name');a=p.parse_args();sys.exit(run(a.name,a.tag,a.case.resolve(),a.output_name))
