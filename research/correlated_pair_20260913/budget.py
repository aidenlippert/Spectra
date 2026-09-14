"""Run one named local experiment with durable timing, logs and a hard timeout."""
import argparse
import datetime
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/correlated_pair_20260913/runs'


def run(name, seconds, command):
    if not name or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for c in name):
        raise ValueError('Safe unique experiment name required')
    if not 1 <= seconds <= 900 or not command: raise ValueError('Bounded command required')
    OUT.mkdir(parents=True, exist_ok=True); path=OUT/(name+'.json'); log=OUT/(name+'.log')
    receipt={'name':name,'command':command,'timeout_seconds':seconds,'status':'starting',
             'started_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    with path.open('x') as f: json.dump(receipt,f,indent=2)
    env=os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMBA_NUM_THREADS'):
        env[key]='1'
    start=time.monotonic(); before=resource.getrusage(resource.RUSAGE_CHILDREN)
    try:
        with log.open('x') as stream:
            proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            receipt['pid']=proc.pid
            try:
                code=proc.wait(timeout=seconds);receipt['status']='passed' if code==0 else 'failed'
            except subprocess.TimeoutExpired:
                receipt['status']='timeout';os.killpg(proc.pid,signal.SIGTERM)
                try: code=proc.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
            receipt['exit_code']=code
    except BaseException as err:
        receipt['status']='runner_error';receipt['error']=repr(err)
        raise
    finally:
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(wall_seconds=time.monotonic()-start,
            child_user_seconds=after.ru_utime-before.ru_utime,child_system_seconds=after.ru_stime-before.ru_stime,
            peak_child_RSS_bytes=after.ru_maxrss*(1 if sys.platform=='darwin' else 1024),
            finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
        path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)
    return 0 if receipt['status']=='passed' else 1


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--seconds',type=int,required=True)
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args();cmd=a.command
    if cmd and cmd[0]=='--':cmd=cmd[1:]
    sys.exit(run(a.name,a.seconds,cmd))
