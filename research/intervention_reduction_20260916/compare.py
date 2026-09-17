"""Matched exact direct-versus-compiled checks, including preparation costs."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from research.intervention_reduction_20260916.kernel import CheckedKernel
from research.intervention_reduction_20260916.robust import check_robust,extend_control_neighborhood
from research.intervention_reduction_20260916.trajectory import check_trajectory


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fixture',type=Path);p.add_argument('proposal',type=Path)
    p.add_argument('output',type=Path);p.add_argument('trajectories',type=Path,nargs='+');a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    data=json.loads(a.fixture.read_text());proposal=json.loads(a.proposal.read_text())
    inputs=[json.loads(path.read_text()) for path in a.trajectories]
    direct=[];times=[]
    for t in inputs:
        start=time.monotonic();budgets=[F(t['horizon_atomic_time'])/400000]*len(proposal['controls'])
        direct.append(check_robust(data,proposal,t,budgets));times.append(time.monotonic()-start)
    start=time.monotonic();kernel=CheckedKernel(data,proposal,inputs[0]['observed_spatial_orbital'])
    prep=time.monotonic()-start;del data,proposal
    compiled=[];matches=[]
    ignored={'exact_replay_seconds','peak_RSS_bytes'}
    for t,old in zip(inputs,direct):
        start=time.monotonic();r=check_trajectory(None,None,t,_kernel=kernel)
        budgets=[F(t['horizon_atomic_time'])/400000]*len(kernel.query_proposal['controls'])
        r=extend_control_neighborhood(kernel.query_data,kernel.query_proposal,t,r,budgets)
        compiled.append(time.monotonic()-start)
        equal={k:v for k,v in old.items() if k not in ignored}=={k:v for k,v in r.items() if k not in ignored}
        matches.append(equal)
        if not equal:
            raise AssertionError('Direct and compiled mathematical receipts differ')
    result={'scope':'same proposals, pulses, exact acceptance and robustness budgets on one host; inherited discovery excluded equally',
        'query_count':len(inputs),'direct_seconds':times,'direct_total_seconds':sum(times),
        'compiled_preparation_seconds':prep,'compiled_query_seconds':compiled,
        'compiled_total_seconds':prep+sum(compiled),'batch_speed_ratio':sum(times)/(prep+sum(compiled)),
        'all_mathematical_fields_identical':all(matches),'comparison_repeated_processes':False,
        'limitation':'single sequential timing, not a hardware variability study or comparison with all conventional solvers'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
