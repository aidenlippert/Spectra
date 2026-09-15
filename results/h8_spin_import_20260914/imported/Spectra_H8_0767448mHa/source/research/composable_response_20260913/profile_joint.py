"""Bounded cold-process profiles; never overwrite scientific certificates."""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import time


def child(case,mode):
    start=time.monotonic()
    from research.composable_response_20260913 import joint,spatial_gap
    data,tail,_,first=joint.load_case(case)
    if mode=='discover':
        gc,gs=joint.propose_gap(data,tail);fc,fs=spatial_gap.propose(data,tail,data['modes']//2-1)
        bounds=joint.second_endpoints(data,tail,first,fc,gc)
        cert={'kind':'composed_joint_sector_response_v1','fixture_sha256':joint.digest(data),
            'first_response_sha256':joint.digest(first),'joint_sector':gc,'first_sector':fc,
            'second_response':joint.outer_program(bounds)}
    else:
        cert=json.loads((joint.OUT/f'{case}_joint_response.json').read_text())
    receipt=joint.check(data,tail,first,cert)
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform!='darwin':peak*=1024
    print(json.dumps({'case':case,'mode':mode,'inside_child_wall_seconds':time.monotonic()-start,
        'peak_RSS_bytes':peak,'receipt':receipt,
        'numerical_packages_loaded':[x for x in ('numpy','scipy','cvxpy','pyscf') if x in sys.modules],
        'inherited_upper_amplitudes_read':False}))


def main():
    root=Path(__file__).resolve().parents[2];out=root/'results/composable_response_20260913'
    rows=[];start=time.monotonic()
    for case in ('h6','h8'):
        for mode in ('discover','replay'):
            t=time.monotonic();command=[sys.executable]+(['-S'] if mode=='replay' else [])
            command+=['-m','research.composable_response_20260913.profile_joint','--child',case,mode]
            run=subprocess.run(command,cwd=root,capture_output=True,text=True,timeout=60)
            if run.returncode:raise RuntimeError(run.stderr)
            row=json.loads(run.stdout);row['cold_process_wall_seconds']=time.monotonic()-t
            if mode=='replay' and row['numerical_packages_loaded']:raise AssertionError('Numerical package in exact replay')
            rows.append(row)
            print(json.dumps({k:v for k,v in row.items() if k!='receipt'}),flush=True)
    result={'profiles':rows,'wall_seconds':time.monotonic()-start,'per_process_timeout_seconds':60,
        'scope':'One measured repetition per case and mode. Discovery includes initial exact acceptance; replay is a separate standard-library-only process. These are lower-component costs from frozen Hamiltonian/tail/scalar-target inputs, not a fresh complete many-body solve.'}
    (out/'profiles.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--child',nargs=2);args=parser.parse_args()
    if args.child:child(*args.child)
    else:main()
