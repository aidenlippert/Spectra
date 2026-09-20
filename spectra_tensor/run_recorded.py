"""Run one construction and record failures plus a process exit receipt."""
from pathlib import Path
from time import perf_counter
import argparse,json,sys,traceback
from .api import solve

def main():
    p=argparse.ArgumentParser();p.add_argument('request',type=Path);p.add_argument('out',type=Path);p.add_argument('--bonds',nargs='+',type=int,default=[48,96,160]);p.add_argument('--sweeps',type=int,default=4);p.add_argument('--seed-steps',type=int,default=12);a=p.parse_args();start=perf_counter();code=1
    try:
        result=solve(json.loads(a.request.read_text()),a.out,bonds=tuple(a.bonds),sweeps=a.sweeps,seed_steps=a.seed_steps)
        code=0 if result['status']=='target_met' else 2
    except BaseException:
        traceback.print_exc()
    finally:
        destination=a.out.parent/(a.out.name+'.process.json')
        with destination.open('x') as f:json.dump(dict(exit_code=code,elapsed_seconds=perf_counter()-start),f)
    sys.exit(code)
if __name__=='__main__':main()
